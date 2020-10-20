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
from mapio.geodict import GeoDict
from impactutils.io.cmd import get_command_output


def trim_ocean(grid2D, mask, all_touched=True, crop=False):
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

    tempdir = tempfile.mkdtemp()

    # Get shapes ready
    if type(mask) == str:
        with fiona.open(mask, 'r') as shapefile:
            bbox = (gdict.xmin, gdict.ymin, gdict.xmax, gdict.ymax)
            hits = list(shapefile.items(bbox=bbox))
            features = [feature[1]["geometry"] for feature in hits]
            # hits = list(shapefile)
            # features = [feature["geometry"] for feature in hits]
    elif type(mask) == list:
        features = mask
    else:
        raise Exception('mask is neither a link to a shapefile or a list of \
                        shapely shapes, cannot proceed')

    if len(features) == 0:
        print('No coastlines in ShakeMap area')
        return grid2D

    tempfilen = os.path.join(tempdir, 'temp.bil')
    tempfile1 = os.path.join(tempdir, 'temp.tif')
    tempfile2 = os.path.join(tempdir, 'temp2.tif')
    GDALGrid.copyFromGrid(grid2D).save(tempfilen)
    cmd = 'gdal_translate -a_srs EPSG:4326 -of GTiff %s %s' % \
        (tempfilen, tempfile1)
    rc, so, se = get_command_output(cmd)

    if rc:
        with rasterio.open(tempfile1, 'r') as src_raster:
            out_image, out_transform = rasterio.mask.mask(src_raster, features,
                                                          all_touched=all_touched,
                                                          crop=crop)
            out_meta = src_raster.meta.copy()
            out_meta.update({"driver": "GTiff",
                             "height": out_image.shape[1],
                             "width": out_image.shape[2],
                             "transform": out_transform})
            with rasterio.open(tempfile2, "w", **out_meta) as dest:
                dest.write(out_image)

        newgrid = GDALGrid.load(tempfile2)

    else:
        print(se)
        raise Exception('ocean trimming failed')

    shutil.rmtree(tempdir)
    return newgrid


def quickcut(filename, gdict, tempname=None, extrasamp=5., method='bilinear',
             precise=True, cleanup=True, verbose=False):
    """
    Use gdal to trim a large global file down quickly so mapio can read it
    efficiently. (Cannot read Shakemap.xml files, must save as .bil filrst)

    Args:
        filename (str): File path to original input file (raster).
        gdict (geodict): Geodictionary to cut around and align with.
        tempname (str): File path to desired location of clipped part of
            filename.
        extrasamp (int): Number of extra cells to cut around each edge of
            geodict to have resampling buffer for future steps.
        method (str): If resampling is necessary, method to use.
        precise (bool): If true, will resample to the gdict as closely as
            possible, if False it will just roughly cut around the area of
            interest without changing resolution
        cleanup (bool): if True, delete tempname after reading it back in
        verbose (bool): if True, prints more details

    Returns: New grid2D layer

    Note: This function uses the subprocess approach because ``gdal.Translate``
        doesn't hang on the command until the file is created which causes
        problems in the next steps.
    """
    if gdict.xmax < gdict.xmin:
        raise Exception('quickcut: your geodict xmax is smaller than xmin')

    try:
        filegdict = GDALGrid.getFileGeoDict(filename)
    except:
        try:
            filegdict = GMTGrid.getFileGeoDict(filename)
        except:
            raise Exception('Cannot get geodict for %s' % filename)

    if tempname is None:
        tempdir = tempfile.mkdtemp()
        tempname = os.path.join(tempdir, 'junk.tif')
        deltemp = True
    else:
        tempdir = None
        deltemp = False

    # if os.path.exists(tempname):
    #     os.remove(tempname)
    #     print('Temporary file already there, removing file')

    filegdict = filegdict[0]

    # Get the right methods for mapio (method) and gdal (method2)
    if method == 'linear':
        method2 = 'bilinear'
    if method == 'nearest':
        method2 = 'near'
    if method == 'bilinear':
        method = 'linear'
        method2 = 'bilinear'
    if method == 'near':
        method = 'nearest'
        method2 = 'near'
    else:
        method2 = method

    if filegdict != gdict:
        # First cut without resampling
        tempgdict = GeoDict.createDictFromBox(
            gdict.xmin, gdict.xmax, gdict.ymin, gdict.ymax,
            filegdict.dx, filegdict.dy, inside=True)

        try:
            egdict = filegdict.getBoundsWithin(tempgdict)

            ulx = egdict.xmin - extrasamp * egdict.dx
            uly = egdict.ymax + extrasamp * egdict.dy
            lrx = egdict.xmax + (extrasamp+1) * egdict.dx
            lry = egdict.ymin - (extrasamp+1) * egdict.dy
            
            cmd = 'gdal_translate -a_srs EPSG:4326 -of GTiff -projwin %1.8f \
            %1.8f %1.8f %1.8f -r %s %s %s' % (ulx, uly, lrx, lry, method2,
                                              filename, tempname)
        except:  
            # When ShakeMap is being loaded, sometimes they won't align right
            # because it's already cut to the area, so just load the whole file
            cmd = 'gdal_translate -a_srs EPSG:4326 -of GTiff -r %s %s %s' % \
                (method2, filename, tempname)
        rc, so, se = get_command_output(cmd)
        if not rc:
            raise Exception(se.decode())
        else:
            if verbose:
                print(so.decode())

        newgrid2d = GDALGrid.load(tempname)
        if precise:
            # Resample to exact geodictionary
            newgrid2d = newgrid2d.interpolate2(gdict, method=method)
        if cleanup:
            os.remove(tempname)

        if deltemp:
            shutil.rmtree(tempdir)

    else:
        ftype = GMTGrid.getFileType(filename)
        if ftype != 'unknown':
            newgrid2d = GMTGrid.load(filename)
        elif filename.endswith('.xml'):
            newgrid2d = ShakeGrid.load(filename)
        else:
            newgrid2d = GDALGrid.load(filename)

    return newgrid2d
