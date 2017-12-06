#!/usr/bin/env python

from mapio.multihaz import MultiHazardGrid
import collections
import json

# NEED TO ADD ERROR CATCHING IN CASE FIELDS ARE MISSING
# TAKES A TON OF MEMORY AND A LONG TIME WITH LARGE FILES


def savelayers(grids, filename):
    """
    Save ground failure layers object as a MultiHazard HDF file, preserving
    metadata structures. Must all have the same geodict.
    Args:
        grids: Ground failure layers object.
        filename (str): Path to where you want to save this file.
    """
    layers = collections.OrderedDict()
    metadata = collections.OrderedDict()
    for key in list(grids.keys()):
        layers[key] = grids[key]['grid'].getData()
        metadata[key] = {
            'description': grids[key]['description'],
            'type': grids[key]['type'],
            'label': grids[key]['label']
        }
    origin = {}
    header = {}
    mgrid = MultiHazardGrid(layers, grids[key]['grid'].getGeoDict(),
                            origin,
                            header,
                            metadata=metadata)
    mgrid.save(filename)


def loadlayers(filename):
    """
    Load a MultiHazard HDF file back in as a ground failure layers object in
    active memory (must have been saved for this purpose).
    Args:
        filename (str): Path to layers file.
    """
    mgrid = MultiHazardGrid.load(filename)
    grids = collections.OrderedDict()
    for key in mgrid.getLayerNames():
        grids[key] = {
            'grid': mgrid.getData()[key],
            'description': mgrid.getMetadata()[key]['description'],
            'type': mgrid.getMetadata()[key]['type'],
            'label': mgrid.getMetadata()[key]['label']
        }

    return grids


def metadata2json(model, filename):
    """
    Save metadata from ground failure model in json file

    Args:
        model (dict): Dictionary of single model output formatted like:

            .. code-block:: python

                {
                    'grid': mapio grid2D object,
                    'label': 'label for colorbar and top line of subtitle',
                    'type': 'output or input to model',
                    'description': 'description for subtitle'
                }
        filename: Name of json file, should end in .json extension
        
    Returns:
        Json file containing metadata information
    """
    metadata = model['description']
    metajson = json.dumps(metadata)
    with open(filename, mode='w') as f3:
        f3.write(metajson)