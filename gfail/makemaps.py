#!/usr/bin/env python

# stdlib imports
import os
import matplotlib as mpl
import tempfile

# from configobj import ConfigObj

# third party imports
import numpy as np
import matplotlib.pyplot as plt
import simplekml
from folium.utilities import mercator_transform

# Make fonts readable and recognizable by illustrator
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = \
    ['Helvetica', 'Arial', 'Bitstream Vera Serif', 'sans-serif']
plt.switch_backend('agg')

# # hex versions:
# DFCOLORS = [
#     '#efefb34D',  # 30% opaque 4D
#     '#e5c72f66',  # 40% opaque 66
#     '#ea720780',  # 50% opaque 80
#     '#c0375c99',  # 60% opaque 99
#     '#5b28b299',  # 60% opaque 99
#     '#1e1e6499'   # 60% opaque 99
# ]
DFCOLORS = [
    [0.94, 0.94, 0.70, 0.7],
    [0.90, 0.78, 0.18, 0.7],
    [0.92, 0.45, 0.03, 0.7],
    [0.75, 0.22, 0.36, 0.7],
    [0.36, 0.16, 0.70, 0.7],
    [0.12, 0.12, 0.39, 0.7]
]

DFBINS = [0.002, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]


def create_kmz(maplayer, outfile, mask=None, levels=None, colorlist=None):
    """
    Create kmz files of models

    Args:
        maplayer (dict): Dictionary of one model result formatted like:

            .. code-block:: python

                {
                    'grid': mapio grid2D object,
                    'label': 'label for colorbar and top line of subtitle',
                    'type': 'output or input to model',
                    'description': 'description for subtitle'
                }
        outfile (str): File extension
        mask (float): make all cells below this value transparent
        levels (array): list of bin edges for each color, must be same length
        colorlist (array): list of colors for each bin, should be length one less than levels

    Returns:
        kmz file
    """
    # Figure out lims
    if levels is None:
        levels = DFBINS
    if colorlist is None:
        colorlist = DFCOLORS

    if len(levels)-1 != len(colorlist):
        raise Exception('len(levels) must be one longer than len(colorlist)')

    # Make place to put temporary files
    temploc = tempfile.TemporaryDirectory()

    # Figure out file names
    name, ext = os.path.splitext(outfile)
    basename = os.path.basename(name)
    if ext != '.kmz':
        ext = '.kmz'
    filename = '%s%s' % (name, ext)
    mapfile = os.path.join(temploc.name, '%s.tiff' % basename)
    legshort = '%s_legend.png' % basename
    legfile = os.path.join(temploc.name, legshort)

    # Make colored geotiff
    out = make_rgba(maplayer['grid'], mask=mask,
                    levels=levels, colorlist=colorlist)
    rgba_img, extent, lmin, lmax, cmap = out
    # Save as a tiff
    plt.imsave(mapfile, rgba_img, vmin=lmin, vmax=lmax, cmap=cmap)

    # Start creating kmz
    L = simplekml.Kml()

    # Set zoom window
    doc = L.document  # have to put lookat in root document directory
    doc.altitudemode = simplekml.AltitudeMode.relativetoground
    boundaries1 = get_zoomextent(maplayer['grid'])
    doc.lookat.latitude = np.mean([boundaries1['ymin'], boundaries1['ymax']])
    doc.lookat.longitude = np.mean([boundaries1['xmax'], boundaries1['xmin']])
    doc.lookat.altitude = 0.
    doc.lookat.range = (boundaries1['ymax']-boundaries1['ymin']) * 111. * 1000.  # dist in m from point
    doc.description = 'USGS near-real-time earthquake-triggered %s model for \
                       event id %s' % (maplayer['description']['parameters']\
                       ['modeltype'], maplayer['description']['event_id'])

    prob = L.newgroundoverlay(name=maplayer['label'])
    prob.icon.href = 'files/%s.tiff' % basename
    prob.latlonbox.north = extent[3]
    prob.latlonbox.south = extent[2]
    prob.latlonbox.east = extent[1]
    prob.latlonbox.west = extent[0]
    L.addfile(mapfile)

    # Add legend and USGS icon as screen overlays
    # Make legend
    make_legend(levels, colorlist, filename=legfile, title=maplayer['label'])
    
    size1 = simplekml.Size(x=0.3, xunits=simplekml.Units.fraction)
    leg = L.newscreenoverlay(name='Legend', size=size1)
    leg.icon.href = 'files/%s' % legshort
    leg.screenxy = simplekml.ScreenXY(x=0.2, y=0.05, xunits=simplekml.Units.fraction,
                                      yunits=simplekml.Units.fraction)
    L.addfile(legfile)

    size2 = simplekml.Size(x=0.15, xunits=simplekml.Units.fraction)
    icon = L.newscreenoverlay(name='USGS', size=size2)
    icon.icon.href = 'files/USGS_ID_white.png'
    icon.screenxy = simplekml.ScreenXY(x=0.8, y=0.95, xunits=simplekml.Units.fraction,
                                       yunits=simplekml.Units.fraction)
    L.addfile(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           os.pardir, 'content', 'USGS_ID_white.png'))

    L.savekmz(filename)
    return filename


def get_zoomextent(grid, propofmax=0.3):
    """
    Get the extent that contains all values with probabilities exceeding
    a threshold in order to determine ideal zoom level for interactive map
    If nothing is above the threshold, uses the full extent

    Args:
        grid: grid2d of model output
        propofmax (float): Proportion of maximum that should be fully included
            within the bounds.

    Returns:
        * boundaries: a dictionary with keys 'xmin', 'xmax', 'ymin', and
         'ymax' that defines the zoomed boundaries in geographic coordinates.

    """
    maximum = np.nanmax(grid.getData())

    xmin, xmax, ymin, ymax = grid.getBounds()
    lons = np.linspace(xmin, xmax, grid.getGeoDict().nx)
    lats = np.linspace(ymax, ymin, grid.getGeoDict().ny)

    if maximum <= 0.:
        boundaries1 = dict(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
        # If nothing is above the threshold, use full extent
        return boundaries1

    threshold = propofmax * maximum

    row, col = np.where(grid.getData() > float(threshold))
    lonmin = lons[col].min()
    lonmax = lons[col].max()
    latmin = lats[row].min()
    latmax = lats[row].max()

    boundaries1 = {}

    if xmin < lonmin:
        boundaries1['xmin'] = lonmin
    else:
        boundaries1['xmin'] = xmin
    if xmax > lonmax:
        boundaries1['xmax'] = lonmax
    else:
        boundaries1['xmax'] = xmax
    if ymin < latmin:
        boundaries1['ymin'] = latmin
    else:
        boundaries1['ymin'] = ymin
    if ymax > latmax:
        boundaries1['ymax'] = latmax
    else:
        boundaries1['ymax'] = ymax

    return boundaries1   


def make_rgba(grid2D, levels, colorlist, mask=None,
              mercator=False):
    """
    Make an rgba (red, green, blue, alpha) grid out of raw data values and
    provide extent and limits needed to save as an image file
    
    Args:
        grid2D: Mapio Grid2D object of result to mape 
        levels (array): list of bin edges for each color, must be same length
        colorlist (array): list of colors for each bin, should be length one
            less than levels
        mask (float): mask all values below this value
        mercator (bool): project to web mercator (needed for leaflet, not
                 for kmz)

    Returns:
        tuple: (rgba_img, extent, lmin, lmax, cmap), where:
            * rgba_img: rgba (red green blue alpha) image
            * extent: list of outside corners of image,
                [minlat, maxlat, minlon, maxlon]
            * lmin: lowest bin edge
            * lmax: highest bin edge
            * cmap: colormap corresponding to image
    """

    data1 = grid2D.getData()
    if mask is not None:
        data1[data1 < mask] = float('nan')
    geodict = grid2D.getGeoDict()
    extent = [
        geodict.xmin - 0.5*geodict.dx,
        geodict.xmax + 0.5*geodict.dx,
        geodict.ymin - 0.5*geodict.dy,
        geodict.ymax + 0.5*geodict.dy,
    ]

    lmin = levels[0]
    lmax = levels[-1]
    data2 = np.clip(data1, lmin, lmax)
    cmap = mpl.colors.ListedColormap(colorlist)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    data2 = np.ma.array(data2, mask=np.isnan(data1))
    rgba_img = cmap(norm(data2))
    if mercator:
        rgba_img = mercator_transform(
            rgba_img, (extent[2], extent[3]), origin='upper')

    return rgba_img, extent, lmin, lmax, cmap


def make_legend(levels, colorlist, filename=None, orientation='horizontal',
                title=None, transparent=False):
    """Make legend file

    Args:

        levels (array): list of bin edges for each color, must be same length
        colorlist (array): list of colors for each bin, should be length one
            less than levels
        filename (str): File extension of legend file
        orientation (str): orientation of colorbar, 'horizontal' or 'vertical'
        title (str): title of legend (usually units)
        transparent (bool): if True, background will be transparent

    Returns:
        figure of legend

    """
    fontsize = 16
    labels = ['< %1.1f%%' % (levels[0] * 100.,)]
    for db in levels:
        if db < 0.01:
            labels.append('%1.1f' % (db * 100,))
        else:
            labels.append('%1.0f' % (db * 100.,))

    if orientation == 'vertical':
        # Flip order to darker on top
        labels = labels[::-1]
        colors1 = colorlist[::-1]
        fig, axes = plt.subplots(len(colors1) + 1, 1,
                                 figsize=(3., len(colors1)-1.7))
        clearind = len(axes)-1
        maxind = 0
    else:
        colors1 = colorlist
        fig, axes = plt.subplots(1, len(colorlist) + 1,
                                 figsize=(len(colorlist) + 1.7, 0.8))
        # DPI = fig.get_dpi()
        # fig.set_size_inches(440/DPI, 83/DPI)
        clearind = 0
        maxind = len(axes)-1

    for i, ax in enumerate(axes):
        ax.set_ylim((0., 1.))
        ax.set_xlim((0., 1.))
        # draw square
        if i == clearind:
            color1 = colors1[0]
            color1[-1] = 0.  # make completely transparent
            if orientation == 'vertical':
                label = labels[i+1]
            else:
                label = labels[0]
        else:
            if orientation == 'vertical':
                label = '%s-%s%%' % (labels[i+1], labels[i])
                color1 = colors1[i]
            else:
                label = '%s-%s%%' % (labels[i], labels[i+1])
                color1 = colors1[i-1]
            color1[-1] = 0.8  # make less transparent
            if i == maxind:
                label = '> %1.0f%%' % (levels[-2]*100.)
        ax.set_facecolor(color1)
        if orientation == 'vertical':
            ax.text(1.1, 0.5, label, fontsize=fontsize,
                    rotation='horizontal', va='center')
        else:
            ax.set_xlabel(label, fontsize=fontsize,
                          rotation='horizontal')

        ax.set_yticks([])
        ax.set_xticks([])
        plt.setp(ax.get_yticklabels(), visible=False)
        plt.setp(ax.get_xticklabels(), visible=False)

    if orientation == 'vertical':
        fig.suptitle(title.title(), weight='bold', fontsize=fontsize+2)
        plt.subplots_adjust(hspace=0.01, right=0.4, top=0.82)
    else:
        fig.suptitle(title.title(), weight='bold', fontsize=fontsize+2)
        # , left=0.01, right=0.99, top=0.99, bottom=0.01)
        plt.subplots_adjust(wspace=0.1, top=0.6)
    # plt.tight_layout()
    if filename is not None:
        fig.savefig(filename, bbox_inches='tight', transparent=transparent)
    else:
        plt.show()


if __name__ == '__main__':
    pass
