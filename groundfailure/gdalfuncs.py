#!/usr/bin/env python
"""
Functions for manipulating geospatial data by building and then calling command
line GDAL functions. GDAL needs to be installed for these to work. Cleans up
temporary files.
"""

import glob
import os

from mapio.gdal import GDALGrid


def getfilesfromfolders(folderwildcard, filewildcard):
    """
    SRTM datasets come in folders, this pulls out the .bil filenames from each
    folder with the full file path

    :param folderwildcard: wild card indicating how the folders are named that
      contain the SRTM .bil files
    :returns:
        filenames: a list of files

    """
    foldernames = glob.glob(folderwildcard)
    filenames = []
    for fold in foldernames:
        tfile = glob.glob('%s/%s' % (fold, filewildcard))
        if len(tfile) > 0:
            filenames += tfile
    return filenames


def srtm2slope(filenames, finaloutfile, fmt='EHdr', cleanup=True):
    """
    Convert tiles of srtm to slope (with intermediate projection to transverse
    mercator projection)
    USAGE srtm2slope(filenames, finaloutfile, fmt='EHdr')

    :param filenames: Single file name as string, Wildcard to call filenames
      [full file path needed (e.g. '/Users/Bob/n*_w*3arc*.bil')] or list of file
      names, should all be in same coordinate system (WGS84 is coordinate system
      for SRTM data from Earth Explorer). Creates intermediate files that can
      be cleaned up or not
    :type filenames: string or list
    :param finaloutfile: Full file path to output file with extension (if no
      path given, will be put in current directory)
    :type finaloutfile: string
    :param fmt: GDAL format to use for outputs, currently only EHdr (.bil) and
      GMT (.grd) are supported, but others might work
    :type fmt: string
    :param cleanup: Whether to delete intermediate files or not
    :type cleanup: Boolean
    :returns:
        a geospatial file (format fmt) named finaloutfile defining slope for the area corresponding to the input filenames
    """

    if '*' in filenames:
        filenames = glob.glob(filenames)
    if fmt == 'EHdr':
        fileext = '.bil'
    elif fmt == 'GMT':
        fileext = '.grd'
    if type(filenames) is str:  # If just one file, don't need to merge
        s_srs, t_srs = warp(filenames, 'temputm58' + fileext, 'EPSG:4326', None,
                            method='bilinear')
    else:
        merge(filenames, 'tempmerge58' + fileext)
        s_srs, t_srs = warp('tempmerge58' + fileext, 'temputm58' + fileext,
                            'EPSG:4326', None, method='bilinear')
    slope('temputm58' + fileext, 'tempslope58' + fileext)
    warp('tempslope58' + fileext, finaloutfile, t_srs, s_srs, method='bilinear')
    # clean up files
    if cleanup:
        delfiles = glob.glob('tempmerge58.*') + glob.glob('temputm58.*') + \
            glob.glob('tempslope58.*')
        for delf in delfiles:
            os.remove(delf)


def merge(filenames, outfile, fmt='EHdr'):
    """Calls gdal_merge to combine tiles (must all be in same coordinate system)
    USAGE merge(filenames, outfile, fmt='EHdr')

    :param filenames: List of full file path to files that will be merged
    :param outfile: Name of output file with extension
    :param fmt: Output format 'EHdr' or 'GMT'
    :returns:
        outfile: a geospatial file of merged tiles
    """

    build = 'gdal_merge.py -of %s -o %s' % (fmt, outfile)
    for filen in filenames:
        build += ' %s' % filen
    #run code
    #subprocess.call(build)
    os.system(build)


def warp(infile, outfile, s_srs, t_srs, method='bilinear', fmt='EHdr'):
    """
    Calls gdalwarp to reproject a raster
    USAGE s_srs, t_srs = warp(infile, outfile, s_srs, t_srs, method='bilinear',
    fmt='EHdr')

    :param infile: Single input file
    :param outfile: Single output file with extension
    :param s_srs: Projection of input file (EPSG or PROJ.4 string), if None,
      uses transverse mercator
    :param t_srs: Projection of output file (EPSG or PROJ.4 string), if None,
      uses transverse mercator
    :param method: method to use in warping, from gdalwarp's options,
      'bilinear', 'nearest' etc., bilinear should be used for slopes to avoid
      weird artifacts
    :param fmt: Format of output, 'EHdr' or 'GMT'
    :returns:
       s_srs: s_srs that was used
       t_srs: t_srs that was used
    """
    if s_srs is None or t_srs is None:
        temp = GDALGrid.getFileGeoDict(infile)
        clat = temp.ymin + (temp.ymax-temp.ymin)/2.0
        clon = temp.xmin + (temp.xmax-temp.xmin)/2.0
    if s_srs is None:
        s_srs = '"+proj=tmerc +lat_0=%s +lon_0=%s +x_0=0 +y_0=0 +units=m +no_defs"' % (clat, clon)
    elif t_srs is None:
        t_srs = '"+proj=tmerc +lat_0=%s +lon_0=%s +x_0=0 +y_0=0 +units=m +no_defs"' % (clat, clon)
    build = 'gdalwarp -overwrite -s_srs %s -t_srs %s -r %s -of %s %s %s' % (s_srs, t_srs, method, fmt, infile, outfile)
    #run code
    os.system(build)
    return s_srs, t_srs


def slope(infile, outfile, v2h=1.0, fmt='EHdr'):
    """
    Calls gdaldem to compute slopes
    USAGE slope(infile, outfile, v2h=1.0, fmt='EHdr')

    :param infile: Single input file
    :type infile: string
    :param outfile: Single output file with extension
    :type outfile: string
    :param v2h: vertical to horizontal scale
    :type v2h: float
    :param fmt: Format of output, 'EHdr' or 'GMT'
    :returns:
        outfile: a geospatial file (format fmt) of slope angle in degrees
    """
    build = 'gdaldem slope %s %s -s %1.0f -of %s' % (infile, outfile, v2h, fmt)
    #run code
    os.system(build)
