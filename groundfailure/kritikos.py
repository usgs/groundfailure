#!/usr/bin/env python

#lib imports
import os.path
import warnings
import urllib.request, urllib.error, urllib.parse
import tempfile
import collections
import math
import numpy as np

#local imports
from mapio.shake import ShakeGrid
from mapio.gdal import GDALGrid
from mapio.geodict import GeoDict


def kritikos_fuzzygamma(shakefile, config, bounds=None):
    """
    Runs kritikos procedure with fuzzy gamma overlay method
    """

    cmodel = config['kritikos_2015']
    gamma = cmodel['gamma_value']

    ## Read in layer files and get data
    layers = cmodel['layers']
    try:
        # Slope
        slope_file = layers['slope']
        # DFF
        dff_file = layers['dff']
        # DFS
        dfs_file = layers['dfs']
    except:
        print('Unable to retrieve grid data.')

    try:
        div = cmodel['divisor']
        # Load in divisors
        MMI_div = div['MMI']
        slope_div = div['slope']
        dff_div = div['dff']
        dfs_div = div['dfs']
        slope_pos_div = div['slope_pos']
    except:
        print('Unable to retrieve divisors.')

    try:
        power = cmodel['power']
        # Load in powers
        MMI_power = power['MMI']
        slope_power = power['slope']
        dff_power = power['dff']
        dfs_power = power['dfs']
        slope_pos_power = power['slope_pos']
    except:
        print('Unable to retrieve powers.')

    # Cut and resample all files
    try:
        shkgdict = ShakeGrid.getFileGeoDict(shakefile, adjust='res')
        slopedict = GDALGrid.getFileGeoDict(slope_file)
        if bounds is not None:  # Make sure bounds are within ShakeMap Grid
            if shkgdict.xmin > bounds['xmin'] or shkgdict.xmax < bounds['xmax'] or shkgdict.ymin > bounds['ymin'] or shkgdict.ymax < bounds['ymax']:
                print('Specified bounds are outside shakemap area, using ShakeMap bounds instead')
                bounds = None
        if bounds is not None:
            tempgdict = GeoDict({'xmin': bounds['xmin'], 'ymin': bounds['ymin'], 'xmax': bounds['xmax'], 'ymax': bounds['ymax'], 'dx': 100., 'dy': 100., 'nx': 100., 'ny': 100.}, adjust='res')
            gdict = slpdict.getBoundsWithin(tempgdict)
        else:  # Get boundaries from shakemap if not specified
            shkgdict = ShakeGrid.getFileGeoDict(shakefile, adjust='res')
            slpdict = GDALGrid.getFileGeoDict(slopefile)
            gdict = slpdict.getBoundsWithin(shkgdict)
    except:
        print('Unable to create base geodict.')

    # Load in data
    try:
        # Load in slope data
        slopegrid = GDALGrid.load(slopefile, samplegeodict=gdict, resample=False)
        slope_data = slopefrid.getData().astype(float)
        # Load in MMI
        shakemap = ShakeGrid.load(shakefile, samplegeodict=gdict, resample=True, method='linear', adjust='res')
        MMI_data = shakemap.getLayer('MMI').getData().astype(float)
        # Load in Dff
        dffgrid = GDALGrid.load(slopefile, samplegeodict=gdict, resample=False)
        dff_data = dffgrid.getData().astype(float)
        # Load in DFS
        dfsgrid = GDALGrid.load(slopefile, samplegeodict=gdict, resample=False)
        dfs_data = dfsgrid.getData().astype(float)
        # Load in Slope Position
        slope_pos_grid = GDALGrid.load(slopefile, samplegeodict=gdict, resample=False)
        slope_pos_data = slop_pos_grid.getData().astype(float)
    except:
        print('Data could not be retrieved.')

    try:
        mmi_class = cmodel['classification']['MMI']
        slope_class = cmodel['classification']['slope']
        dff_class = cmodel['classification']['dff']
        dfs_class = cmodel['classification']['dfs']
        slope_pos_class = cmodel['classification']['slope_pos']
    except:
        print('Could not recover classifications from config.')

    try:
        mmi_classes = mmi_class.split(',')
        for i in mmi_classes:
            if i.find('-') != -1:
                j = i.split('-')
                if MMI_data in range(int(j[0]), int(j[1])):
                    MMI_data = int(j[0])
            else:
                MMI_data = int(i)
    except:
        print('Could not categorize MMI values')

    try:
        slope_classes = slope_class.split(',')
        k = 1
        for i in mmi_classes:
            if i.find('-') != -1:
                j = i.split('-')
                if slope_data in range(int(j[0]), int(j[1])):
                    slope_data = k
                    k += 1
            else:
                slope_data = 11
    except:
        print('Could not recategorize Slope Values.')

    try:
        dff_classes = dff_class.split(',')
        k = 1
        for i in dff_classes:
            if i.find('-') != -1:
                j = i.split('-')
                if dff_data in range(int(j[0]), int(j[1])):
                    dff_data = k
                    k += 1
            else:
                dff_data = 7
    except:
        print('Could not recategorize DFF values.')

    try:
        dfs_classes = dfs_class.split(',')
        k = 1
        for i in dfs_classes:
            if i.find('-') != -1:
                j = i.split('-')
                if dfs_data in range(int(j[0]), int(j[1])):
                    dfs_data = k
                    k += 1
            else:
                dfs_data = 6
    except:
        print('Could not recategorize DFS values.')

    try:
        slope_pos_classes = slope_pos_class.split(',')
        k = 1
        for i in slope_poss_classes:
            if slope_pos_data == i:
                slope_pos_data = k
                k += 1
    except:
        print('Could not recategorize slope position values.')

    try:
        layers = []
        # Calculate layers
        slope = 1/(1+np.exp(slope_data/slope_div, slope_power))
        MMI = 1/(1+np.exp(MMI_data/MMI_div, MMI_power))
        dff = 1/(1+np.exp(dff_data/dff_div, dff_power))
        dfs = 1/(1+np.exp(dfs_data/dfs_div, dfs_power))
        slope_pos = 1/(1+np.exp(slop_pos_data/slop_pos_div, slope_pos_power))
        # Add to layers list (to be used in further calculations)
        layers.append(slope)
        layers.append(MMI)
        layers.append(dff)
        layers.append(dfs)
        layers.append(slope_pos)
    except:
        print('Layer calculations failed.')

    try:
        # Calculate final model
        for l in layers.items():


        mu_x = (pi_operator(mu_i, 1, n))^(1-gamma) * (1 - pi_operator(1-mu_i, 1, n))^(gamma)
        # pi_operator is a multiplication summation












