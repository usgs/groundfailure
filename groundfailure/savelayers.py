#!/usr/bin/env python

from mapio.multihaz import MultiHazardGrid
import collections

#NEED TO ADD ERROR CATCHING IN CASE FIELDS ARE MISSING


def savelayers(grids, filename):
    """
    Save ground failure layers object as a MultiHazard HDF file, preserving metadata structures
    Must all have the same geodict
    :param grids: ground failure layers object
    :param filename: Full file path to where you want to save this file
    """
    layers = collections.OrderedDict()
    metadata = collections.OrderedDict()
    for key in grids.keys():
        layers[key] = grids[key]['grid']
        metadata[key] = {'description': grids[key]['description'], 'type': grids[key]['type'], 'label': grids[key]['label']}
    origin = {}
    header = {}
    mgrid = MultiHazardGrid(layers, grids[key]['grid'].getGeoDict(), origin, header, metadata=metadata)
    mgrid.save(filename)


def loadlayers(filename):
    """
    Load a MultiHazard HDF file back in as a ground failure layers object in active memory (must have been saved for this purpose)
    """
    mgrid = MultiHazardGrid.load(filename)
    grids = collections.OrderedDict()
    for key in mgrid.getLayerNames():
        grids[key]['grid'] = mgrid.getData()[key]
        grids[key]['description'] = mgrid.getMetadata()[key]['description']
        grids[key]['type'] = mgrid.getMetadata()[key]['type']
        grids[key]['label'] = mgrid.getMetadata()[key]['label']

    return grids
