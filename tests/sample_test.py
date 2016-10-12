#!/usr/bin/env python

# Tests of sample.py

# stdlib imports
import numpy as np
from groundfailure import sample
from collections import OrderedDict
from mapio import Grid2D


def _test_sample():
    shapes = []
    poly1 = [(34.25, 34.5),
             (34.3, 34.4),
             (34.4, 34.4),
             (34.45, 34.5),
             (34.5, 34.7),
             (34.4, 34.7),
             (34.3, 34.65),
             (34.25, 34.6)]
    poly2 = [(34.75, 34.65),
             (34.75, 34.35),
             (34.8, 34.35),
             (34.85, 34.3),
             (34.95, 34.3),
             (34.95, 34.6)]
    shp1 = {'id': 0,
            'type': 'Feature',
            'properties': OrderedDict([('value', 1)]),
            'geometry': {'type': 'Polygon', 'coordinates': poly1}}
    shp2 = {'id': 1,
            'type': 'Feature',
            'properties': OrderedDict([('value', 2)]),
            'geometry': {'type': 'Polygon', 'coordinates': poly2}}
    shapes.append(shp1)
    shapes.append(shp2)
    allverts = np.vstack((np.array(poly1), np.array(poly2)))
    xmin = allverts[:, 0].min()
    xmax = allverts[:, 0].max()
    ymin = allverts[:, 1].min()
    ymax = allverts[:, 1].max()
    bounds = (xmin, ymin, xmax, ymax)
    dx = 5500.0
    YesTestPoints, YesTrainPoints, NoTestPoints, NoTrainPoints, xvar, yvar, pshapes, proj = sample.sampleFromShapes(shapes, bounds, dx=dx, Nsamp=50, testPercent=0.5)
    sample.plotPoints(shapes, YesTestPoints, YesTrainPoints, NoTestPoints, NoTrainPoints, 'output.png')

    #make up a geology vector data set
    geopoly1 = [(34.0, 34.6),
                (34.0, 34.0),
                (34.1, 34.0),
                (34.15, 34.1),
                (34.3, 34.25),
                (34.35, 34.4),
                (34.25, 34.5),
                (34.25, 34.6)]
    geopoly2 = [(34.25, 34.6),
                (34.25, 34.5),
                (34.35, 34.4),
                (34.6, 34.4),
                (34.6, 34.55),
                (34.6, 34.6)]
    geopoly3 = [(34.1, 34.0),
                (34.85, 34.0),
                (34.6, 34.4),
                (34.35, 34.4),
                (34.3, 34.25),
                (34.15, 34.1)]
    geopoly4 = [(34.85, 34.0),
                (35.3, 34.0),
                (35.3, 34.6),
                (34.6, 34.6),
                (34.6, 34.4),
                (34.85, 34.15)]
    value = 4
    idx = 0
    shapes = []
    for poly in [geopoly1, geopoly2, geopoly3, geopoly4]:
        shp = {'id': idx,
               'type': 'Feature',
               'properties': OrderedDict([('value', value)]),
               'geometry': {'type': 'Polygon', 'coordinates': poly}}
        shapes.append(shp)
        value += 1
        idx += 1

    #make up a grid data set
    geodict = {'xmin': 33.5, 'xmax': 35.5, 'ymin': 33.5, 'ymax': 35.0, 'dx': 0.5, 'dy': 0.5, 'ny': 4, 'nx': 5}
    data = np.arange(0, 20).reshape((4, 5))
    grid = Grid2D(data, geodict)

    yestestgeovalues = sample.sampleShapes(shapes, YesTestPoints, 'value')
    yestestgridvalues = sample.sampleFromGrid(grid, YesTestPoints)

    notestgeovalues = sample.sampleShapes(shapes, NoTestPoints, 'value')
    notestgridvalues = sample.sampleFromGrid(grid, YesTestPoints)

    TestFrame = pd.DataFrame()
    yescov = np.ones_like(yestestgeovalues[:, 0])
    nocov = np.ones_like(notestgeovalues[:, 0])
    yeslat = yestestgeovalues[:, 1]
    yeslon = yestestgeovalues[:, 0]
    nolat = notestgeovalues[:, 1]
    nolon = notestgeovalues[:, 0]
    TestFrame['Latitude'] = np.vstack((yeslat, nolat))
    TestFrame['Longitude'] = np.vstack((yeslat, nolat))
    TestFrame['Coverage'] = np.vstack((yescov, nocov))
    TestFrame['Geology'] = np.vstack((yestestgeovalues, notestgeovalues))
