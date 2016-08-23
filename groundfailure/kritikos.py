#!/usr/bin/env python

"""
Kritikos based landslide statistic_models
"""

"""
These methods requires:
    1) a comprehenive and accurate landslide inventory
    2) multiple strong motion recordings or detailed isoseismals
    3) accurate geologic maps
    4) extensive geotechnical data
    5) accurate digital elevation models (DEMs)
Utilizes semi data-driven memberships
"""


"""
Note that none of the functions used the ‘Hedge’ option. We combined these using the fuzzy Gamma Overlay with a gamma value of 0.9.

Density of landslides calculated from a frequency ratio

Frequency Ration = (N_Li / N_Ci) / (N_L / N_A)                      (3)

where   N_Li is the number of landslide pixels in the factor i
        N_Ci is the total number of pixels in the factor i
        N_L is total number of landslide pixels in the study area
        N_A is the total number of pixels of the study area
"""

"""
Operators = fuzzyAND, fuzzyOR, fuzzyPRODUCT, fuzzySUM, fuzzyGAMMA

fuzzyGAMMA operator ----
mu_x = (pi_operator(mu_i, 1, n))^(1-gamma) * (1 - pi_operator(1-mu_i, 1, n))^(gamma)
# pi_operator is a multiplication summation
"""

############# Config file
'''
[statistic_models]
    [[kritikos_2014]]
        gamma_value = 0.9
        [[[layers]]]
            slope = Users/kbiegel/Documents/GroundFailure/inputs/Northridge/Northridge_SLP_WGS84_1arcsec.bil
            dff = 
            dfs = 
            slope_pos =
        [[[divisor]]]
            MMI = 7.5
            slope = 4.875
            dff = 2.375
            dfs = 3.25
            slope_pos = 2.325
        [[[power]]]
            MMI = -14
            slope = -2.65
            dff = 5.375
            dfs = 5.5
            slope_pos = -4.375
        [[[classification]]]
            MMI = 5,6,7,8,9
            slope = 0-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49, 50+  # Reclassify as 1,2,3,etc.
            dff = 0-4, 5-9, 10-19, 20-29, 30-39, 40-49, 50+  # Reclassify as 1,2,3,etc.
            dfs = 0-0.49, 0.5-0.99, 1.0-1.49, 1.5-1.99, 2.0-2.49, 2.5+  # Reclassify as 1,2,3,etc.
            slope_pos = 'Flat', 'Valley', 'Mid-Slope', 'Ridge'  # Reclassify as 1,2,3,etc.
'''

config = {'statistic_models': {'kritikos_2014': {'gamma_value': 0.9, 'layers': {'slope': 'Users/kbiegel/Documents/GroundFailure/inputs/Northridge/Northridge_SLP_WGS84_1arcsec.bil', 'dff': , 'dfs': , 'slope_pos': }, 'divisor': {'MMI': 7.5, 'slope': 4.875, 'dff': 2.375, 'dfs': 3.25, 'slope_pos': 2.325}, 'power': {'MMI': -14, 'slope': -2.65, 'dff': 5.375, 'dfs': 5.5, 'slope_pos': -4.375}, 'classification': {'MMI': [5,6,7,8,9], 'slope': [range(0,4), range(5,9), range(10,14), range(15,19), range(20,24), range(25,29), range(30,34), range(35,39), range(40,44), range(45,49), 50], 'dff': [range(0,4), range(5,9), range(10,19), range(20,29), range(30,39), range(40,49), 50], 'dfs': [[0,0.49],[0.5,0.99], [1.0,1.49], [1.5,1.99], [2.0,2.49], 2.5], 'slope': ['Flat', 'Valley', 'Mid-Slope', 'Ridge']}}}}

######## Config File

#stdlib imports
import os.path
import warnings
import urllib.request, urllib.error, urllib.parse
import tempfile
import collections

#turn off all warnings...
warnings.filterwarnings('ignore')

#local imports
from mapio.shake import ShakeGrid
from mapio.gdal import GDALGrid
from mapio.geodict import GeoDict

#third party imports
import numpy as np


def kritikos_fuzzygamma(shakefile, config, bounds=None):
    """
    Runs kritikos procedure with fuzzy gamma
    """

    cmodel = config['statistic_models']['kritikos_2014']
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
        # Slope Position
        slope_pos_file = layers['slope_pos']
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

            [[[classification]]]
            MMI = 5,6,7,8,9
            slope = 0-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49, 50+  # Reclassify as 1,2,3,etc.
            dff = 0-4, 5-9, 10-19, 20-29, 30-39, 40-49, 50+  # Reclassify as 1,2,3,etc.
            dfs = 0-0.49, 0.5-0.99, 1.0-1.49, 1.5-1.99, 2.0-2.49, 2.5+  # Reclassify as 1,2,3,etc.
            slope_pos = 'Flat', 'Valley', 'Mid-Slope', 'Ridge'  # Reclassify as 1,2,3,etc.

    try:
        ranges = {}
        for keys in cmodel['classification']:




    try:
        if MMI_data is in range(5,6):
            MMI_data = 5
        elif MMI_data is in range(6,7):
            MMI_data = 6
        elif MMI_data is in range(7,8):
            MMI_data = 7
        elif MMI_data is in range(8,9):
            MMI_data = 8
        else:
            MMI_data = 9
    except:
        print('Could not categorize MMI values')

    try:
        if slope_data is in range(0,4):
            slope_data = 1
        elif slope_data is in range(5,9):
            slope_data = 2
        elif slope_data is in range(10,14):
            slope_data = 3
        elif slope_data is in range(15,19):
            slope_data = 4
        elif slope_data is in range(20,24):
            slope_data = 5
        elif slope_data is in range(25,29):
            slope_data = 6
        elif slope_data is in range(30,34):
            slope_data = 7
        elif slope_data is in range(35,39):
            slope_data = 8
        elif slope_data is in range(40,44):
            slope_data = 9
        elif slope_data is in range(45,49):
            slope_data = 10
        else:
            slope_data = 11
    except:
        print('Could not recategorize Slope Values.')

    try:
        if dff_data is in range(0,4):
            dff_data = 1
        elif dff_data is in range(5,9):
            dff_data = 2
        elif dff_data is in range(10,19):
            dff_data = 3
        elif dff_data is in range(20,29):
            dff_data = 4
        elif dff_data is in range(30,39):
            dff_data = 5
        elif dff_data is in range(40,49):
            dff_data = 6
        else:
            dff_data = 7
    except:
        print('Could not recategorize DFF values.')

    try:
        if int(dfs_data/10.) is in range(0,5):
            dfs_data = 1
        elif int(dfs_data/10) is in range(5,10):
            dfs_data = 2
        elif int(dfs_data/10) is in range(10,15):
            dfs_data = 3
        elif int(dfs_data/10) is in range(15,20):
            dfs_data = 4
        elif int(dfs_data/10) is in range(20,25):
            dfs_data = 5
        else:
            dfs_data = 6
    except:
        print('Could not recategorize DFS values.')

    try:
        if slope_pos_data is 'Flat':
            slope_pos_data = 1
        elif slope_pos_data is 'Valley':
            slope_pos_data = 2
        elif slope_pos_data is 'Mid-Slope':
            slope_pos_data = 3
        elif slope_pos_data is 'Ridge':
            slope_pos_data = 4
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












