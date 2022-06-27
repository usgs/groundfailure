#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spatial functions used by all models
"""

import os
import tempfile
import fiona
import shutil
import rasterio
import rasterio.mask

from mapio.gdal import GDALGrid
from mapio.shake import ShakeGrid
from mapio.gmt import GMTGrid
from mapio.grid2d import Grid2D
from mapio.geodict import GeoDict, geodict_from_affine
from mapio.writer import write
from impactutils.io.cmd import get_command_output
import numpy as np


def split_grid(grid):
    # split grid that sits on 180 meridian
    xmin, xmax, ymin, ymax = (
        grid._geodict.xmin,
        grid._geodict.xmax,
        grid._geodict.ymin,
        grid._geodict.ymax,
    )
    dx, dy = grid._geodict.dx, grid._geodict.dy

    toprow, topcol = grid._geodict.getRowCol(ymax, 180.0)
    leftdata = grid._data[:, 0:topcol]
    leftny, leftnx = leftdata.shape
    leftdict = {
        "xmin": xmin,
        "xmax": 180.0,
        "ymin": ymin,
        "ymax": ymax,
        "dx": dx,
        "dy": dy,
        "nx": leftnx,
        "ny": leftny,
    }
    leftgeodict = GeoDict(leftdict)
    leftgrid = Grid2D(data=leftdata, geodict=leftgeodict)

    rightdata = grid._data[:, topcol:]
    rightny, rightnx = rightdata.shape
    rightdict = {
        "xmin": -180.0,
        "xmax": xmax,
        "ymin": ymin,
        "ymax": ymax,
        "dx": dx,
        "dy": dy,
        "nx": rightnx,
        "ny": rightny,
    }
    rightgeodict = GeoDict(rightdict)
    rightgrid = Grid2D(data=rightdata, geodict=rightgeodict)

    return (leftgrid, rightgrid)


def join_grids(leftgrid, rightgrid):
    leftdict = leftgrid.getGeoDict()
    rightdict = rightgrid.getGeoDict()
    newdata = np.concatenate((leftgrid._data, rightgrid._data), axis=1)
    ny, nx = newdata.shape
    newdict = {
        "xmin": leftdict.xmin,
        "xmax": rightdict.xmax,
        "ymin": leftdict.ymin,
        "ymax": leftdict.ymax,
        "dx": leftdict.dx,
        "dy": leftdict.dy,
        "nx": nx,
        "ny": ny,
    }
    geodict = GeoDict(newdict)
    newgrid = Grid2D(data=newdata, geodict=geodict)
    return newgrid


def trim_ocean2(grid2D, mask, all_touched=True, crop=False):
    """Use the mask (a shapefile) to trim offshore areas

    Args:
        grid2D: MapIO grid2D object of results that need trimming
        mask: list of shapely polygon features already loaded in or string of
            file extension of shapefile to use for clipping
        all_touched (bool): if True, won't mask cells that touch any part of
            polygon edge
        crop (bool): crop boundaries of raster to new masked area

    Returns:
        grid2D file with ocean masked
    """
    gdict = grid2D.getGeoDict()

    # Get shapes ready
    if type(mask) == str:
        with fiona.open(mask, "r") as shapefile:
            if gdict.xmin < gdict.xmax:
                bbox = (gdict.xmin, gdict.ymin, gdict.xmax, gdict.ymax)
                hits = list(shapefile.items(bbox=bbox))
                features = [feature[1]["geometry"] for feature in hits]
            else:
                features = []
                bbox1 = (gdict.xmin, gdict.ymin, 180.0, gdict.ymax)
                hits = list(shapefile.items(bbox=bbox1))
                tfeatures = [feature[1]["geometry"] for feature in hits]
                features += tfeatures
                bbox2 = (-180, gdict.ymin, gdict.xmax, gdict.ymax)
                hits = list(shapefile.items(bbox=bbox2))
                tfeatures = [feature[1]["geometry"] for feature in hits]
                features += tfeatures
    elif type(mask) == list:
        features = mask
    else:
        raise Exception(
            "mask is neither a link to a shapefile or a list of \
                        shapely shapes, cannot proceed"
        )

    if len(features) == 0:
        print("No coastlines in ShakeMap area")
        return grid2D

    # now make a dataset from our Grid2D object.
    # I think the only way to do this is to write out a file
    # and then open it. Stupid.
    try:
        tempdir = tempfile.mkdtemp()
        tmpfile = os.path.join(tempdir, "tmp.tif")
        # TODO: Create a temp dir, write a file in it
        if gdict.xmin < gdict.xmax:
            write(grid2D, tmpfile, "tiff")
            dataset = rasterio.open(tmpfile, "r")
            out_image, out_transform = rasterio.mask.mask(
                dataset, features, all_touched=all_touched, nodata=np.nan, crop=crop
            )
            dataset.close()
            out_image = np.squeeze(out_image)
            ny, nx = out_image.shape
            geodict = geodict_from_affine(out_transform, ny, nx)
            newgrid = Grid2D(data=out_image, geodict=geodict)
        else:
            leftgrid, rightgrid = split_grid(grid2D)
            write(leftgrid, tmpfile, "tiff")
            dataset = rasterio.open(tmpfile, "r")
            out_image, out_transform = rasterio.mask.mask(
                dataset, features, all_touched=all_touched, nodata=np.nan, crop=crop
            )
            dataset.close()
            out_image = np.squeeze(out_image)
            ny, nx = out_image.shape
            geodict = geodict_from_affine(out_transform, ny, nx)
            newleftgrid = Grid2D(data=out_image, geodict=geodict)

            write(rightgrid, tmpfile, "tiff")
            dataset = rasterio.open(tmpfile, "r")
            out_image, out_transform = rasterio.mask.mask(
                dataset, features, all_touched=all_touched, nodata=np.nan, crop=crop
            )
            dataset.close()
            out_image = np.squeeze(out_image)
            ny, nx = out_image.shape
            geodict = geodict_from_affine(out_transform, ny, nx)
            newrightgrid = Grid2D(data=out_image, geodict=geodict)

            newgrid = join_grids(newleftgrid, newrightgrid)

        return newgrid
    except Exception as e:
        msg = f'Failure to trim input grid: "{str(e)}".'
        raise Exception(msg)
    finally:
        # TODO: clean up
        shutil.rmtree(tempdir)
