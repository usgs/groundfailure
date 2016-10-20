#!/usr/bin/env python

#stdlib imports
import os
import fiona
import pyproj
from functools import partial
from shapely.ops import transform
from shapely.geometry import shape, Point, MultiPoint
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
from scipy import interpolate
from sklearn.metrics import roc_curve, roc_auc_score, auc
import copy

#local imports
from groundfailure.sample import pointsFromShapes
from mapio.gdal import GDALGrid

# Make fonts readable and recognizable by illustrator
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = ['Arial', 'Bitstream Vera Serif', 'sans-serif']


def modelSummary(models, titles=None, outputtype=None, plottype='hist', bounds=None, bins=10, semilogy=False, thresh=0., excludenon=False, showplots=True, saveplots=False, filepath=None):
    """
    Function for creating a summary histogram of a model output

    :param model: Grid2D object of model results or list of Grid2D objects of multiple models
    :param titles: List of titles to use for each model, in same order as model. If none, will use key of each model (may be non-unique)
    :param outputtype: Type of model output, just used for label (e.g.'probability', 'coverage', 'index')
    :param plottype: 'hist' or 'cumul' for histogram or cumulative plot
    :param bounds: Bounding box to include in summary as dictionary e.g. {'xmin': -119.2, 'xmax': -118., 'ymin': 34., 'ymax': 34.7}. If None, will use entire area in model
    :param bins: bins to use for histogram and pie chart, if a single integer, will create that number of bins using min and max of data, if a numpy array, will use those as bin edges
    :param semilogy: = uses log scale instead of linear on y axis if True
    :param thresh: threshold for a nonoccurrence, default is zero but for models that never have nonzero values, can set to what you decide is insignificant
    :param excludenon: If True, will exclude cells that are <=thresh from plot results
    :param showplots: if True, will display the plots
    :param saveplots: if True, will save the plots
    :param filepath: Filepath for saved plots, if None, will save in current directory. Files are named with test name and time stamp
    :returns: Grid2D object of difference between inventory and model (inventory - model)
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)

    means = []
    vallist = []
    if type(models) != list:
        models = [models]
    if titles is None:
        titles = [list(mod.keys()) for mod in models]
    for k, mod in enumerate(models):
        grid = mod.getData()
        if type(bins) is int:
            bins = np.linspace(grid.min(), grid.max(), bins)
        # Trim model, if needed
        if bounds is not None and len(bounds) == 4:
            model1 = copy.deepcopy(mod)  # to avoid changing original file
            model1 = model1.cut(bounds['xmin'], bounds['xmax'], bounds['ymin'], bounds['ymax'], align=True)
        else:
            model1 = mod

        grid = model1.getData()
        allvals = grid[~np.isnan(grid)]
        means.append(np.mean(allvals))
        total = len(allvals)
        totalnonz = len(allvals[allvals > float(thresh)])
        if excludenon is True:
            totalf = totalnonz
            allvals = allvals[allvals > float(thresh)]
        else:
            totalf = total

        # if plottype == 'pie':
        #     hist, bin_edges = np.histogram(allvals, bins=bins)
        #     labels = ['%0.1f-%0.1f' % (bin_edges[i], bin_edges[i+1]) for i in range(len(bin_edges)-1)]
        #     output = ax.pie(hist/float(totalf)*100, labels=None, autopct='%1.1f%%', shadow=False, startangle=90)
        #     plt.legend(output[0], labels, loc="best")
        #     plt.tight_layout()
        #     #import pdb; pdb.set_trace()
        #     plt.axis('equal')

        if plottype == 'cumul':
            values, base = np.histogram(allvals, bins=200)
            cumulative = np.cumsum(values)
            ax.plot(base[:-1], cumulative/float(totalf), label='%s - %1.1e' % (titles[k], means[k]))
            #sd = np.sort(allvals)
            #ax.step(sd, np.arange(sd.size)/float(totalf))
            ax.set_ylabel('Proportion of cells (cumulative)')
        else:
            vallist.append(allvals)
            if k == len(models)-1:
                n, bins, rects = ax.hist(tuple(vallist), bins=bins, normed=True)
                # for rect in rects:
                #     height = rect.get_height()
                #     plt.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%0.1f%%' % (height/float(totalf),), ha='center', va='bottom')
                ax.set_ylabel('Proportion of cells')
        ax.set_xlabel(outputtype)

        # Put a legend below current axis
        ax.legend(loc=9, bbox_to_anchor=(0.5, -0.1), fancybox=True, ncol=2)

    if semilogy is True:
        ax.set_yscale('log')

    # if excludenon is True:
    #     plt.suptitle('%d cells > %0.1f out of %d total cells (%0.2f%%)\nTotals plotted exclude cells <= %0.1f' % (totalnonz, thresh, total, totalnonz/float(total), thresh))
    # else:
    #     plt.suptitle('%d cells > %0.1f out of %d total cells (%0.2f%%)' % (totalnonz, thresh, total, totalnonz/float(total)))
    #plt.tight_layout()
    #plt.subplots_adjust(top=0.95)
    plt.subplots_adjust(bottom=0.25)

    if saveplots is True:
        if filepath is None:
            filepath = os.getcwd()
        import datetime
        time1 = datetime.datetime.utcnow().strftime('%d%b%Y_%H%M')
        fig.savefig(os.path.join(filepath, 'Summary_%s.pdf' % (time1,)))

    if showplots is True:
        plt.show()
    else:
        plt.close(fig)

    return ax


def computeCoverage_accurate(gdict, inventory, numdiv=10.):
    """
    VERY SLOW!!
    Slow but more accurate method to produce grid of area actually affected by landsliding in each cell defined by geodict

    :param gdict: geodict, likely taken from model to compare inventory against
    :param inventory: full file path to shapefile of inventory, must be in geographic coordinates, WGS84
    :type inventory: string
    :param numdiv: Approximate amount to subdivide each cell of geodict by to compute areas (higher number slower but more accurate)

    :returns: Grid2D object reporting areal coverage of landsliding inside each cell defined by geodict
    """

    f = fiona.collection(inventory, 'r')
    shapes = list(f)
    bxmin, bymin, bxmax, bymax = f.bounds

    lons = np.linspace(gdict.xmin, gdict.xmax, gdict.nx)
    lats = np.linspace(gdict.ymax, gdict.ymin, gdict.ny)
    llons, llats = np.meshgrid(lons, lats)

    spacing = np.round(np.abs(((lats[1]-lats[0])*111.12*1000.)/numdiv))  # in meters
    yespoints, nopoints, xvar, yvar, pshapes, proj = pointsFromShapes(shapes, bounds=(gdict.xmin, gdict.ymin, gdict.xmax, gdict.ymax), dx=spacing)

    # Loop over lat lon pairs that are within boundaries of yes and no points
    ptlonmax = (np.max((yespoints[:, 0].max(), nopoints[:, 0].max())))
    ptlonmin = (np.max((yespoints[:, 0].min(), nopoints[:, 0].min())))
    ptlatmax = (np.max((yespoints[:, 1].max(), nopoints[:, 1].max())))
    ptlatmin = (np.max((yespoints[:, 1].min(), nopoints[:, 1].min())))

    subllons = llons[(llons >= ptlonmin) & (llons <= ptlonmax) & (llats >= ptlatmin) & (llats <= ptlatmax)]
    subllats = llats[(llons >= ptlonmin) & (llons <= ptlonmax) & (llats >= ptlatmin) & (llats <= ptlatmax)]

    import time
    # Contains points method
    t1 = time.clock()
    dx = gdict.dx
    area = np.zeros(np.shape(llons))
    numpts = area.copy()
    numyes = area.copy()
    for lat1, lon1 in zip(subllats, subllons):
        # Find ratio of yes points to no points
        bbPath = mplPath.Path(np.array([[lon1-0.5*dx, lat1-0.5*dx], [lon1-0.5*dx, lat1+0.5*dx], [lon1+0.5*dx, lat1+0.5*dx], [lon1+0.5*dx, lat1-0.5*dx]]))
        yesin = sum(bbPath.contains_points(yespoints))  # sum([(yes0 > lon1-0.5*dx) & (yes0 <= lon1+0.5*dx) & (yes1 > lat1-0.5*dx) & (yes1 <= lat1+0.5*dx)])
        noin = sum(bbPath.contains_points(nopoints))  # sum([(no0 > lon1-0.5*dx) & (no0 <= lon1+0.5*dx) & (no1 > lat1-0.5*dx) & (no1 <= lat1+0.5*dx)])
        total = yesin + noin
        if total == 0.:
            continue
        # get indices
        row = np.where(lats == lat1)
        col = np.where(lons == lon1)
        # Store total number of points in matrix
        numpts[row, col] = total
        # Store area
        numyes[row, col] = yesin
    t2 = time.clock()
    print(('Time elapsed %0.2f seconds' % (t2-t1)))

    # Correct for incompletely sampled squared (all unsampled points would be no points)
    numpts[numpts < (numpts[numpts != 0].mean() - numpts[numpts != 0].std())] = np.median(numpts[numpts != 0])  # Will change zeros to nonzeros, but yeses will be 0 in those cells so it doesn't matter
    area = numyes/numpts

    inventorygrid = GDALGrid(area, gdict)

    return inventorygrid


def computeCoverage(gdict, inventory, numdiv=30., method='nearest'):
    """Fast method to produce grid of area actually affected by landsliding in each cell defined by geodict

    :param gdict: geodict, likely taken from model to compare inventory against
    :param inventory: full file path to shapefile of inventory, must be in geographic coordinates, WGS84
    :type inventory: string
    :param numdiv: Approximate amount to subdivide each cell of geodict by to compute areas (higher number slower but more accurate)
    :return inventorygrid: Grid2D object reporting proportional area of landsliding inside each cell defined by geodict
    :param method: method for resampling when projecting back to geographic coordinates, nearest recommended but not perfect. Cubic not recommended.

    :returns: Grid2D object reporting approximate areal coverage of input inventory corresponding to geodict
    """

    f = fiona.collection(inventory, 'r')
    shapes = list(f.items(bbox=(gdict.xmin, gdict.ymin, gdict.xmax, gdict.ymax)))  # get only shapes that are intersect area of interest from gdict
    bxmin, bymin, bxmax, bymax = f.bounds

    lons = np.linspace(gdict.xmin, gdict.xmax, gdict.nx)
    lats = np.linspace(gdict.ymax, gdict.ymin, gdict.ny)

    # # Convert Shapes and input grid to Transverse mercator
    bclat = bymin + (bymax-bymin)/2.0
    bclon = bxmin + (bxmax-bxmin)/2.0
    original = pyproj.Proj(f.crs)
    destination = pyproj.Proj(projparams='+proj=tmerc +datum=WGS84 +lat_0=%0.5f +lon_0=%0.5F +units=meters +x_0=0 +y_0=0' % (bclat, bclon))
    project = partial(pyproj.transform, original, destination)
    unproject = partial(pyproj.transform, destination, original)
    pshapes = []
    for tshape in shapes:
        gshape = shape(tshape[1]['geometry'])
        pshape = transform(project, gshape)
        pshapes.append(pshape)

    # Get corners for sampling mesh from gdict
    urcnr = Point(gdict.xmax, gdict.ymax)
    llcnr = Point(gdict.xmin, gdict.ymin)
    urcnr1 = transform(project, urcnr)
    llcnr1 = transform(project, llcnr)

    # Figure out what max extent of shapes is
    sxmax = 0.
    sxmin = 1.e12
    symax = 0.
    symin = 1.e12

    for pshape in pshapes:
        x, y = np.array(pshape.exterior.coords.xy)
        if x.max() > sxmax:
            sxmax = x.max()
        if x.min() < sxmin:
            sxmin = x.min()
        if y.max() > symax:
            symax = y.max()
        if y.min() < symin:
            symin = y.min()

    # Create mesh in projected space that would oversample by numdiv
    numsampx = (np.array(urcnr.coords.xy[0]) - np.array(llcnr.coords.xy[0]))/gdict.dx * numdiv
    numsampy = (np.array(urcnr.coords.xy[1]) - np.array(llcnr.coords.xy[1]))/gdict.dy * numdiv

    xdiff = sxmax - sxmin
    ydiff = symax - symin
    xvec = np.linspace(np.array(llcnr1.coords.xy[0]), np.array(urcnr1.coords.xy[0]), numsampx)
    xvec = xvec[(xvec > sxmin-0.15*xdiff) & (xvec < sxmax+0.15*xdiff)]
    yvec = np.linspace(np.array(llcnr1.coords.xy[1]), np.array(urcnr1.coords.xy[1]), numsampy)
    yvec = yvec[(yvec > symin-0.15*ydiff) & (yvec < symax+0.15*ydiff)]

    values = np.zeros((len(yvec), len(xvec)))

    # Basically rasterize. There could be a better way to do this using rasterio, but hard to maintain the right grids
    shapeidx = 0
    yespoints = []
    for pshape in pshapes:
        if not shapeidx % 2000:
            print('Searching polygon %i of %i' % (shapeidx, len(pshapes)))
        shapeidx += 1
        pxmin, pymin, pxmax, pymax = pshape.bounds
        try:
            leftcol = np.where((pxmin - xvec) >= 0)[0].argmax()
            rightcol = np.where((xvec - pxmax) >= 0)[0][0]
            bottomrow = np.where((pymin - yvec) >= 0)[0].argmax()
            toprow = np.where((yvec - pymax) >= 0)[0][0]
        except:  # Shape out of bounds, skip
            continue
        xp = xvec[leftcol:rightcol+1]
        yp = yvec[bottomrow:toprow+1]
        xmesh, ymesh = np.meshgrid(xp, yp)
        xy = list(zip(xmesh.flatten(), ymesh.flatten()))
        for point in xy:
            ix = np.where(xvec == point[0])[0][0]
            iy = np.where(yvec == point[1])[0][0]
            if pshape.contains(Point(point)):
                yespoints.append(point)
                values[iy, ix] = 1

    # Block mean to downsample
    areaval = block_mean(values, numdiv)
    xvecnew = line_mean(xvec, numdiv)
    yvecnew = line_mean(yvec, numdiv)

    xm, ym = np.meshgrid(xvecnew, yvecnew)
    nzpts = MultiPoint(list(zip(*[xm[areaval > 0], ym[areaval > 0]])))
    nzvals = areaval[areaval > 0]

    # Project nonzero points back
    nzpts_proj = transform(unproject, nzpts)

    # Find nearest coordinate for each on original grid

    import time
    t1 = time.clock()
    if method == 'nearest':  # This is more efficient than using nearest with griddata
        area = np.zeros((len(lats), len(lons)))
        latidx = 0
        for pt, val in zip(nzpts_proj, nzvals):
            if not latidx % 1000:
                print('Searching coord %i of %i' % (latidx, len(nzpts_proj)))
            latidx += 1
            lon, lat = np.array(pt.coords.xy)
            row = (np.abs(lats - lat)).argmin()
            col = (np.abs(lons - lon)).argmin()
            area[row, col] = val
        t2 = time.clock()
    else:
        latidx = 0
        PTS = np.squeeze(np.array([pt.coords.xy for pt in nzpts_proj]))
        llons, llats = np.meshgrid(lons, lats)
        numpts = np.shape(llons)[0]*np.shape(llons)[1]
        xdatpts = np.reshape(llons, numpts)
        ydatpts = np.reshape(llats, numpts)
        area = interpolate.griddata(PTS, nzvals, list(zip(xdatpts, ydatpts)), method=method.lower(), fill_value=0.)
        area = np.reshape(area, np.shape(llons))
        t2 = time.clock()
    print(('Time elapsed %0.2f seconds' % (t2-t1)))

    inventorygrid = GDALGrid(area, gdict)

    return inventorygrid


def statsCoverage(modelgrid, inventorygrid, bins=None, showplots=True, saveplots=False, filepath=None):
    """TO DO - FIND MORE TESTS THAT ARE EASY TO COMPARE WITH EACH OTHER
    Compute stats and make comparison plots specific to models that output areal coverage like Godt et al 2008

    :param modelgrid: Grid2D object of model results
    :param inventorygrid: Grid2D object of areal coverage of inventory computed on same grid as modelgrid using, for example, computeCoverage
    :param bins: bin edges to use for various binning and threshold statistical calculations. if None bins = np.linspace(0, np.max((inv.max(), model.max())), 11)
    :param showplots: if True, will display the plots
    :param saveplots: if True, will save the plots
    :param filepath: Filepath for saved plots, if None, will save in current directory. Files are named with test name and time stamp
    :returns:
        * results: dictionary of results of tests.
                      {'Compare coverage': dictionary,
                        'RMS': float,
                        'RMS_nonzero': float,
                        'Percent in bin': dictionary,
                        'Direct Comparison': dictionary]}
        * invminusmod: Grid2D object of difference between inventory and model (inventory - model)

    """

    inv = inventorygrid.getData()
    model = modelgrid.getData()
    invminusmod = GDALGrid(inv-model, inventorygrid.getGeoDict())

    plt.ioff()

    # Make statistical comparisons
    results = {}
    if bins is None:
        bins = np.linspace(0, np.max((inv.max(), model.max())), 11)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist((model, inv), bins=bins)
    ax.set_yscale('log')
    ax.legend(('Model', 'Inventory'))
    ax.set_ylabel('Total # of cells')
    ax.set_xlabel('Predicted coverage bin')
    results['Direct Comparison'] = {'bins': bins, 'model': model, 'inventory': inv}

    # Statistical comparison
    perc = []
    areatot = []
    for i, bin in enumerate(bins[:-1]):
        idx = (model > bin) & (model <= bins[i+1])
        areain = inv[idx]
        totalin = sum((areain > bin) & (areain <= bins[i+1]))
        perc.append(float(totalin)/len(idx)*100)
        areatot.append(np.mean(areain))
    areatot = np.nan_to_num(areatot)

    binvec = []
    for i in range(len(bins[:-1])):
        binvec.append(bins[i]+(bins[i+1]-bins[i])/2)
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(binvec, perc, 'o-')
    ax1.set_ylabel('% in right bin')
    ax1.set_xlabel('Predicted coverage bin')
    ax1.set_title('Actual coverage in predicted bin')
    results['Percent in bin'] = {'bincenters': binvec, 'perc': perc}

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(binvec, areatot, 'o-')
    ax2.plot([0., 1.], [0., 1.], '--', color='gray')
    ax2.set_ylabel('Actual coverage')
    ax2.set_xlabel('Predicted coverage')
    ax2.set_title('Actual vs. Predicted coverage')
    results['Compare coverage'] = {'bincenters': binvec, 'areatot': areatot}

    if showplots is True:
        plt.show()
    if saveplots is True:
        if filepath is None:
            filepath = os.getcwd()
        import datetime
        time1 = datetime.datetime.utcnow().strftime('%d%b%Y_%H%M')
        fig.savefig(os.path.join(filepath, 'Direct_compare_%s.pdf' % (time1,)))
        fig1.savefig(os.path.join(filepath, 'Percent_in_bin_%s.pdf' % (time1,)))
        fig2.savefig(os.path.join(filepath, 'Compare_coverage_%s.pdf' % (time1,)))

    # RMS
    idx = np.union1d(model.nonzero()[0], inv.nonzero()[0])
    results['RMS'] = np.sqrt((model - inv)**2).mean()
    print(('RMS: %0.3f' % results['RMS']))
    results['RMS_nonzero'] = np.sqrt((model[idx] - inv[idx])**2).mean()
    print(('RMS_nonzero: %0.3f' % results['RMS_nonzero']))

    return invminusmod, results


def stats(modelgrid, inventory, dx=100., Nsamp=None, method='nearest', extent='inventory', bins=None, runtests=True, showplots=True, saveplots=False, filepath=None):
    """Run through suite of tests for models that output probability or index that varies between 0 and 1

    :param modelgrid: Grid2D object of model results
    :param inventory: full file path to shapefile of inventory, must be in geographic coordinates, WGS84
    :type inventory: string
    :param dx: Approximate sample spacing in meters, overwritten if Nsamp is defined
    :type dx: float
    :param Nsamp: Total number of samples desired - will choose optimal dx to get slightly more than this number of samples delete samples outside bounds and randomly delete others until sample number is exactly Nsamp
    :type Nsamp: integer
    :param method: method used for interp2d when transforming sampled model values back from projected coordinates to geographic coordinates - 'nearest', 'linear', or 'cubic'
    :type method: string
    :param extent: extent to include in sampling - 'inventory' or 'model' or custom bounds as tuple (xmin, ymin, xmax, ymax) - in lats and lons
    :param bins: bin edges to use for various binning and threshold statistical calculations. if None bins = [0, 0.2, 0.4, 0.6, 0.8, 1.]
    :param runtests: if True, will run various statistical tests, if False will just output sampled values
    :param showplots: if True, will disply the plots
    :param saveplots: if True, will save the plots
    :param filepath: Filepath for saved plots, if None, will save in current directory. Files are named with test name and time stamp
    :returns:
        * yespoints: Nx2 array of geographic coordinates of positive sample locations
        * nopoints: Nx2 array of geographic coordinates of negative sample locations
        * modelvalyes: N model output values corresponding to yespoints
        * modelvalno: N model output values corresponding to nopoints
        * results: dictionary of results of statistical tests. Will be empty if runtests=False
        {'Occ_nonocc': dict,
        'SRC': dict,
        'ROC': dict,
        'AUC_ROC': float,
        'Log_loss': float,
        'GFC': dict,
        'Pred_vs_Obs': dict,
        'Brier': float,
        'Brier_no': float,
        'Brier_yes': float}
    """
    plt.close('all')
    f = fiona.collection(inventory, 'r')
    shapes = list(f)
    bxmin, bymin, bxmax, bymax = f.bounds
    gdict = modelgrid.getGeoDict()

    if extent == 'model':
        extent = gdict.xmin, gdict.ymin, gdict.xmax, gdict.ymax
    elif extent == 'inventory':
        extent = bxmin, bymin, bxmax, bymax
    #yespoints, junk, nopoints, junk2, xvar, yvar, pshapes, proj = sampleFromShapes(shapes, extent, dx=dx, Nsamp=Nsamp, testPercent=100.)

    yespoints, nopoints, xvar, yvar, pshapes, proj = pointsFromShapes(shapes, extent, dx=dx, Nsamp=Nsamp)
    yesptx = [pt[0] for pt in yespoints]
    yespty = [pt[1] for pt in yespoints]
    noptx = [pt[0] for pt in nopoints]
    nopty = [pt[1] for pt in nopoints]
    #import pdb; pdb.set_trace()
    # Get values of model at those points
    lons = np.linspace(gdict.xmin, gdict.xmax, gdict.nx)
    lats = np.linspace(gdict.ymax, gdict.ymin, gdict.ny)
    if method.lower() == 'nearest':
        modelvalyes = []
        modelvalno = []
        for XX, YY in zip(yesptx, yespty):
            row = (np.abs(lats - YY)).argmin()
            col = (np.abs(lons - XX)).argmin()
            modelvalyes.append(modelgrid.getData()[row, col])
        for XX, YY in zip(noptx, nopty):
            row = (np.abs(lats - YY)).argmin()
            col = (np.abs(lons - XX)).argmin()
            modelvalno.append(modelgrid.getData()[row, col])
    else:
        func = interpolate.interp2d(lons, lats, modelgrid.getData(), kind=method.lower())
        modelvalyes = np.array([float(func(XX, YY)) for XX, YY in zip(yesptx, yespty)])
        modelvalno = np.array([float(func(XX, YY)) for XX, YY in zip(noptx, nopty)])

    modelvalyes = np.nan_to_num(np.array(modelvalyes))  # replace nan with zeros
    modelvalno = np.nan_to_num(np.array(modelvalno))  # replace nan with zeros

    # Now run the desired tests and make the desired plots
    results = {}

    if runtests is True:
        # Brier score
        N = len(yespoints) + len(nopoints)
        yessum = np.sum([(val-1)**2 for val in modelvalyes])
        nosum = np.sum([(val)**2 for val in modelvalno])
        results['Brier_yes'] = yessum/len(modelvalyes)
        results['Brier_no'] = nosum/len(modelvalno)
        results['Brier'] = (yessum + nosum)/N
        print(('Brier scores: overall %0.3f\nBrier_yes score: %0.3f\nBrier_no score %0.3f' % (results['Brier'], results['Brier_yes'], results['Brier_no'])))

        # Logarithmic score
        tempno = np.array(modelvalno).copy()
        tempyes = np.array(modelvalyes).copy()
        tempno[tempno == 0] = 1.e-15
        tempyes[tempyes == 0] = 1.e-15
        results['Log_loss'] = -(np.sum(np.log(tempyes)) + np.sum(np.log(1.-tempno)))/N
        print(('Log loss score: %0.3f' % (results['Log_loss'],)))

        if bins is None:
            bins = [0, 0.2, 0.4, 0.6, 0.8, 1.]
        binvec = []
        observed = []
        percyes = []
        percno = []
        overall_tot = len(modelvalyes) + len(modelvalno)
        for i in range(len(bins[:-1])):
            binvec.append(bins[i]+(bins[i+1]-bins[i])/2)
            yestot = np.sum([(modelvalyes > bins[i]) & (modelvalyes < bins[i+1])])
            notot = np.sum([(modelvalno > bins[i]) & (modelvalno < bins[i+1])])
            if notot+yestot != 0:
                observed.append(float(yestot)/(yestot+notot))
            else:
                observed.append('nan')
            percyes.append((yestot/float(overall_tot))*100.)
            percno.append((notot/float(overall_tot))*100.)

        plt.ioff()

        # Predicted vs. Observed ratios
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(binvec, observed, '-o')
        ax.plot([0]+binvec, [0]+binvec, '--', color='gray')
        ax.set_xlabel('Expected ratio')
        ax.set_ylabel('Observed ratio')
        ax.set_xlim([bins[0], bins[-1]])
        ax.set_title('Predicted vs. Observed')
        results['Pred_vs_Obs'] = {'binvec': binvec, 'observed': observed}

        # Ground failure occurrence/nonoccurrence
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        wid = (bins[1]-bins[0])/2.5
        rects1 = ax1.bar(np.array(bins[:-1]), percyes, width=wid)
        rects2 = ax1.bar(np.array(bins[:-1])+wid, percno, width=wid, color='r')
        ax1.set_xlabel('Predicted susceptibility range')
        ax1.set_ylabel('% of samples')
        ax1.legend((rects1[0], rects2[0]), ('Occurrence', 'Nonoccurrence'))
        ax1.set_title('Occurrence vs. Nonoccurrence')
        results['Occ_nonocc'] = {'bins': bins, 'percyes': percyes, 'percno': percno}

        # Ground failure capture for various thresholds
        gfc = []
        for val in bins:
            gfc.append(np.sum([modelvalyes > val])/float(len(yespoints)))
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        ax2.plot(bins, gfc, 'o-')
        ax2.set_xlabel('Threshold')
        ax2.set_ylabel(r'%GFC')
        ax2.set_title('Ground Failure Capture')
        results['GFC'] = {'thresholds': bins, 'gfc': gfc}

        # ROC curves
        fpr, tpr, thresholds = roc_curve(np.concatenate((np.ones(len(yespoints)), np.zeros(len(nopoints)))), np.concatenate((modelvalyes, modelvalno)))
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        ax3.plot(fpr, tpr)
        ax3.set_xlabel('False positive rate')
        ax3.set_ylabel('True positive rate')
        ax3.set_xlim([0, 1.])
        ax3.set_ylim([0, 1.])
        ax3.plot(fpr, fpr, '--', color='gray')
        ax3.set_title('ROC curve')
        results['ROC'] = {'thresholds': bins, 'gfc': gfc}

        results['AUC_ROC'] = roc_auc_score(np.concatenate((np.ones(len(yespoints)), np.zeros(len(nopoints)))), np.concatenate((modelvalyes, modelvalno)))
        print(('AUC_ROC: %0.3f' % (results['AUC_ROC'],)))
        ax3.text(0.8, 0.2, 'AUC: %0.3f' % results['AUC_ROC'])

        # Success rate curves
        sucbin = np.linspace(0, 1., 100)
        prop = []
        realvals = np.concatenate((np.ones(len(yespoints)), np.zeros(len(nopoints))))
        predvals = np.concatenate((modelvalyes, modelvalno))
        indx = np.argsort(predvals)
        predvals = predvals[indx]
        realvals = realvals[indx]
        for val in sucbin:
            prop.append(np.sum(realvals[predvals < val])/len(yespoints))
        fig4 = plt.figure()
        ax4 = fig4.add_subplot(111)
        ax4.plot(sucbin, prop)
        ax4.set_xlabel('Success Rate Curve')
        ax4.set_ylabel('Proportion of actual occurrences')
        ax4.set_title('Proportion of Study Area')
        AUC = auc(sucbin, prop)
        print(('AUC_SRC: %0.3f' % AUC))
        ax4.text(0.8, 0.2, 'AUC: %0.3f' % AUC)
        ax4.set_xlim([0, 1.])
        ax4.set_ylim([0, 1.])
        results['SRC'] = {'xvals': sucbin, 'proportion': prop, 'auc': AUC}

        if showplots is True:
            plt.show()
        if saveplots is True:
            if filepath is None:
                filepath = os.getcwd()
            import datetime
            time1 = datetime.datetime.utcnow().strftime('%d%b%Y_%H%M')
            fig.savefig(os.path.join(filepath, 'Pred_vs_obs_%s.pdf' % (time1,)))
            fig1.savefig(os.path.join(filepath, 'Occ_nonocc_%s.pdf' % (time1,)))
            fig2.savefig(os.path.join(filepath, 'GFC_%s.pdf' % (time1,)))
            fig3.savefig(os.path.join(filepath, 'ROC_%s.pdf' % (time1,)))
            fig4.savefig(os.path.join(filepath, 'SRC_%s.pdf' % (time1,)))

    return yespoints, nopoints, modelvalyes, modelvalno, results


def block_mean(ar, fact):
    """
    Block mean for downsampling 2d array
    From here http://stackoverflow.com/questions/18666014/downsample-array-in-python
    """
    from scipy import ndimage
    sx, sy = ar.shape
    newx = int(np.floor(sx/float(fact)))
    newy = int(np.floor(sy/float(fact)))
    vec = np.reshape(np.arange(newx*newy), (newx, newy))
    regions = vec.repeat(fact, 0).repeat(fact, 1)
    # Patch on edges (just expand existing cells on the edges a little)
    xdiff = ar.shape[0] - regions.shape[0]
    if xdiff != 0:
        row = regions[-1, :]
        regions = np.row_stack((regions, np.repeat(np.reshape(row, (1, len(row))), xdiff, axis=0)))
    ydiff = ar.shape[1] - regions.shape[1]
    if ydiff != 0:
        col = regions[:, -1]
        regions = np.column_stack((regions, np.repeat(np.reshape(col, (len(col), 1)), ydiff, axis=1)))
    res = ndimage.mean(ar, labels=regions, index=np.arange(regions.max() + 1))
    res = res.reshape(newx, newy)
    return res


def line_mean(ar, fact):
    """
    Same as block_mean but for 1d array
    """
    from scipy import ndimage
    sx = len(ar)
    newx = int(np.floor(sx/float(fact)))
    vec = np.arange(newx)
    regions = vec.repeat(fact)
    # Patch on edges (just expand existing cells on the edges a little)
    xdiff = ar.shape[0] - regions.shape[0]
    if xdiff != 0:
        row = regions[-1]
        regions = np.concatenate((regions, np.repeat(row, xdiff)))
    res = ndimage.mean(ar, labels=regions, index=np.arange(regions.max() + 1))
    return res
