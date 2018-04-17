#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module contains a set of functions for sampling from inventory
shapefiles for model development or testing. The codes are borrowed
heavily from <http://github.com/mhearne-usgs/lsprocess>
"""

# stdlib imports
from functools import partial

import numpy as np
import pyproj
from shapely.ops import transform
from shapely.geometry import Polygon, shape, MultiPoint, Point

from mapio.geodict import GeoDict
import matplotlib.path as mplPath
from rasterio.transform import Affine
import rasterio


def getProjectedShapes(shapes, xmin, xmax, ymin, ymax):
    """
    Take a sequence of geographic shapes and project them to a bounds-centered
    orthographic projection.

    Args:
        shapes: Sequence of shapes, as read in by fiona.collection().
        xmin: Eastern boundary of all shapes.
        xmax: Western boundary of all shapes.
        ymin: Southern boundary of all shapes.
        ymax: Northern boundary of all shapes.

    Returns:
       Tuple of
           - Input sequence of shapes, projected to orthographic
           - PyProj projection object used to transform input shapes
    """
    latmiddle = ymin + (ymax - ymin) / 2.0
    lonmiddle = xmin + (xmax - xmin) / 2.0
    projstr = ('+proj=ortho +datum=WGS84 +lat_0=%.4f +lon_0=%.4f '
               '+x_0=0.0 +y_0=0.0' % (latmiddle, lonmiddle))
    proj = pyproj.Proj(projparams=projstr)
    project = partial(
        pyproj.transform,
        pyproj.Proj(proj='latlong', datum='WGS84'),
        proj)

    pshapes = []
    for tshape in shapes:
        if tshape['geometry']['type'] == 'Polygon':
            pshapegeo = shape(tshape['geometry'])
        else:
            pshapegeo = shape(tshape['geometry'])
        pshape = transform(project, pshapegeo)
        pshapes.append(pshape)  # assuming here that these are simple polygons

    return (pshapes, proj)


def getYesPoints(pshapes, proj, dx, nmax, touch_center=True):
    """
    Collect x/y coordinates of all points within hazard coverage polygons at
    desired resolution.

    Args:
        pshapes: Sequence of orthographically projected shapes.
        proj: PyProj projection object used to transform input shapes.
        dx: Float resolution of grid at which to sample points, must be round
            number.
        nmax: Threshold maximum number of points in total data mesh.
        touch_center: Boolean indicating whether presence of polygon in each
            grid cell is enough to turn that into a yes pixel.
            ** This doc needs a better explanation! ** And I don't know what
            does so I'm not the one to fix it.

    Returns:
        A tuple of
          - numpy 2-D array of X/Y coordinates inside hazard polygons.
          - number of rows of resulting mesh
          - number of columns of resulting mesh
          - numpy array of x coordinate centers of columns
          - numpy array of y coordinate centers of rows
          - 1D array of indices where yes pixels are located (use np.unravel_index to unpack to 2D array)

    """

    mxmin = 9e10
    mxmax = -9e10
    mymin = 9e10
    mymax = -9e10
    for pshape in pshapes:
        pxmin, pymin, pxmax, pymax = pshape.bounds
        if pxmin < mxmin:
            mxmin = pxmin
        if pxmax > mxmax:
            mxmax = pxmax
        if pymin < mymin:
            mymin = pymin
        if pymax > mymax:
            mymax = pymax

    if not touch_center:
        geodict = GeoDict.createDictFromBox(mxmin, mxmax, mymin, mymax, dx, dx)
        img = rasterizeShapes(pshapes, geodict)
        # now get the numpy array of x/y coordinates where covgrid == 1
        idx = np.where(img == 1)[0]
        x, y = np.unravel_index(idx, (geodict.ny, geodict.nx))
        yespoints = list(zip(x.flatten(), y.flatten()))
        nrows = geodict.ny
        ncols = geodict.nx
        xvar = np.arange(geodict.xmin, geodict.xmax + geodict.dx, geodict.dx)
        yvar = np.arange(geodict.ymin, geodict.ymax + geodict.dy, geodict.dy)
    else:
        xvar = np.arange(mxmin, mxmax + dx, dx)
        yvar = np.arange(mymin, mymax + dx, dx)
        ncols = len(xvar)
        nrows = len(yvar)
        if nmax is not None:
            if ncols * nrows > nmax:
                aspect = ncols / nrows
                ncols = np.sqrt(nmax * aspect)
                nrows = nmax / ncols
                ncols = int(ncols)
                nrows = int(nrows)
                # re-calculate dx here...
                tdx = (mxmax - mxmin) / ncols
                tdy = (mymax - mymin) / nrows
                dx = np.max(tdx, tdy)
                xvar = np.arange(mxmin, mxmax + dx, dx)
                yvar = np.arange(mymin, mymax + dx, dx)

        # Get the "yes" points to sample from
        yespoints = []
        idx = []
        shapeidx = 0
        if pshapes[0].type == 'Polygon':
            # loop over shapes, projecting each one, then get the sample points
            for pshape in pshapes:
                if not shapeidx % 1000:
                    print('Searching polygon %i of %i' %
                          (shapeidx, len(pshapes)))
                shapeidx += 1
                pxmin, pymin, pxmax, pymax = pshape.bounds
                leftcol = np.where((pxmin - xvar) >= 0)[0].argmax()
                rightcol = np.where((xvar - pxmax) >= 0)[0][0]
                bottomrow = np.where((pymin - yvar) >= 0)[0].argmax()
                toprow = np.where((yvar - pymax) >= 0)[0][0]
                xp = np.arange(xvar[leftcol], xvar[rightcol] + dx, dx)
                yp = np.arange(yvar[bottomrow], yvar[toprow] + dx, dx)
                xmesh, ymesh = np.meshgrid(xp, yp)
                xy = list(zip(xmesh.flatten(), ymesh.flatten()))
                for point in xy:
                    ix = np.where(xvar == point[0])[0][0]
                    iy = np.where(yvar == point[1])[0][0]
                    if pshape.contains(Point(point)):
                        yespoints.append(point)
                        idx.append(np.ravel_multi_index(
                            (iy, ix), (nrows, ncols), mode='raise', order='C'))
        else:
            yespoints = []
            for pshape in pshapes:
                yespoints.append(pshape.coords[0])

    return (np.array(yespoints), nrows, ncols, xvar, yvar, idx)


def pointsFromShapes(shapes, bounds, dx=10.0, nmax=None, Nsamp=None,
                     touch_center=True):
    """
    Get yes/no points from shapefile input - same as sampleFromShapes but
    without class balance or separation of test and train, only samples in
    box enclosing the polygons

    Args:
        shapes: Sequence of projected shapes.
        bounds: Tuple of xmin, ymin, xmax, ymax, in lat/lon coordinates, only
            will accept points from within these bounds.
        dx: Resolution of sampling in X and Y (meters), must be a round number
            of meters.
        nmax: If not None, maximum allowed number of mesh points in X and Y
            together (nrows*ncols).  Overrides dx.
        Nsamp: If not None, maximum number of total samples, keeps proportion
            of yes's and no's the same.
        touch_center: Boolean indicating whether presence of polygon in each
            grid cell is enough to turn that into a yes pixel.
            ** This explanation needs to be improved because it seems to
            describe the opposite of the usual meaning of this variable!**

    Returns:
        Tuple of
          - sequence of coordinates in lat/lon for: YesPoints, NoPoints
          - numpy array of mesh column centers
          - numpy array of mesh row centers
          - PyProj object defining orthographic projection of xy points

    """
    xmin, ymin, xmax, ymax = bounds
    shptype = shapes[0]['geometry']['type']
    if shptype not in ['Polygon']:
        raise Exception('Only polygon data types supported!')

    # Get the shapes projected into an orthographic projection centered on
    # the data
    pshapes, proj = getProjectedShapes(shapes, xmin, xmax, ymin, ymax)

    # Get the projected bounds
    project = partial(
        pyproj.transform,
        pyproj.Proj(proj='latlong', datum='WGS84'),
        proj)
    bbPoly = Polygon(((xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin)))
    bbPolyproj = transform(project, bbPoly)

    if Nsamp is not None:
        # Recompute dx, underestimate by dividing by 1.5 so later trimming
        # doesn't reduce desired total
        projbounds = bbPolyproj.bounds
        dx = np.round(np.sqrt(
            ((projbounds[2] - projbounds[0]) *
             (projbounds[3] - projbounds[1])) / (Nsamp)) / 1.5)

    # get the "yes" sample points
    yespoints, nrows, ncols, xvar, yvar, yesidx = getYesPoints(
        pshapes, proj, dx, nmax=nmax, touch_center=touch_center)

    # sampleNo but with taking all of the points instead of just some of them
    # randomly flattened array of all indices in mesh
    allidx = np.arange(0, len(xvar) * len(yvar))
    noidx = np.setxor1d(allidx, yesidx)  # allidx - avoididx
    rowidx, colidx = np.unravel_index(noidx, (len(yvar), len(xvar)))
    nopoints = []
    for row, col in zip(rowidx, colidx):
        xp = xvar[col]
        yp = yvar[row]
        nopoints.append((xp, yp))
    nopoints = np.array(nopoints)

    # Only accept points inside the bounds
    bbPath = mplPath.Path(
        (list(zip(*np.array(bbPolyproj.exterior.coords.xy)))))
    yespoints = yespoints[bbPath.contains_points(yespoints)]
    nopoints = nopoints[bbPath.contains_points(nopoints)]
    totalpoints = (len(nopoints) + len(yespoints))

    if Nsamp is not None and totalpoints > Nsamp:
        ratioyes = float(len(yespoints)) / totalpoints
        keepy = int(np.round(ratioyes * Nsamp))
        indy = np.random.randint(0, len(yespoints), size=keepy)
        indn = np.random.randint(0, len(nopoints), size=Nsamp - keepy)
        yespoints = yespoints[indy, :]
        nopoints = nopoints[indn, :]

    elif totalpoints < Nsamp:
        print(('Only collected %1.0f points out of desired %1.0f points '
               'due to bound restrictions' % (totalpoints, Nsamp)))

    # project all of the point data sets back to lat/lon
    yespoints = projectBack(yespoints, proj)
    nopoints = projectBack(nopoints, proj)

    return (yespoints, nopoints, xvar, yvar, pshapes, proj)


def projectBack(points, proj):
    """
    Project a 2D array of XY points from orthographic projection to decimal
    degrees.

    Args:
        points: 2D numpy array of XY points in orthographic projection.
        proj: PyProj object defining projection.

    Returns:
        2D numpy array of Lon/Lat coordinates.

    """

    mpoints = MultiPoint(points)
    project = partial(
        pyproj.transform,
        proj,
        pyproj.Proj(proj='latlong', datum='WGS84'))
    gmpoints = transform(project, mpoints)
    coords = []
    for point in gmpoints.geoms:
        x, y = point.coords[0]
        coords.append((x, y))
    coords = np.array(coords)
    return coords


def rasterizeShapes(pshapes, geodict, all_touched=True):
    """
    Rasterizing a shape

    Args:
        pshapes: Sequence of orthographically projected shapes.
        goedict: Mapio geodictionary.
        all_touched: Turn pixel "on" if shape touches pixel, otherwise turn it
            on if the center of the pixel is contained within the shape. Note
            that the footprint of the shape is inflated and the amount of
            inflation depends on the pixel size if all_touched=True.

    Returns:
        Rasterio grid.
    """
    outshape = (geodict.ny, geodict.nx)
    txmin = geodict.xmin - geodict.dx / 2.0
    tymax = geodict.ymax + geodict.dy / 2.0
    transform = Affine.from_gdal(
        txmin, geodict.dx, 0.0, tymax, 0.0, -geodict.dy)
    imgs = rasterio.features.rasterize(
        pshapes, out_shape=outshape, fill=0.0, transform=transform,
        all_touched=all_touched, default_value=1.0)
    return imgs
