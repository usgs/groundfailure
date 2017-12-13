#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
from configobj import ConfigObj

# third party
from gfail.conf import correct_config_filepaths
import gfail.logisticmodel as LM
from gfail import stats
import matplotlib
matplotlib.use('Agg')

# where is this script?
homedir = os.path.dirname(os.path.abspath(__file__))
upone = os.path.join(homedir, os.pardir)
datadir = os.path.abspath(os.path.join(homedir, 'data'))


def test_stats_models():
    conf_file = os.path.join(upone, 'defaultconfigfiles', 'models',
                             'zhu_2015.ini')
    conf = ConfigObj(conf_file)
    data_path = os.path.join(datadir, 'loma_prieta', 'model_inputs')
    conf = correct_config_filepaths(data_path, conf)
    shakefile = os.path.join(datadir, 'loma_prieta', 'grid.xml')
    lm = LM.LogisticModel(shakefile, conf, saveinputs=True)
    #maplayers1 = lm.calculate()
    conf_file = os.path.join(upone, 'defaultconfigfiles', 'models',
                             'zhu_2017_coastal.ini')
    conf = ConfigObj(conf_file)
    data_path = os.path.join(datadir, 'loma_prieta', 'model_inputs')
    conf = correct_config_filepaths(data_path, conf)
    shakefile = os.path.join(datadir, 'loma_prieta', 'grid.xml')
    lm = LM.LogisticModel(shakefile, conf, saveinputs=True)
    maplayers2 = lm.calculate()
    # Change shakemap name so that it doesn't stomp on the other
    maplayers2['model']['description']['shakemap'] = '19891018000415_ver2'

    #model_list = [maplayers1, maplayers2]
    #test_dict1 = assess_models.concatenateModels(model_list)
    #test_dict2 = assess_models.concatenateModels(model_list, astitle='model')

    # I think this test is freezing on travis
#    tmp = assess_models.modelSummary(test_dict2, showplots=False,
#                                     summary_figure=False,
#                                     individual_plots=False)
#    np.testing.assert_allclose(tmp[0][0], 0.025677016713957716)
#    np.testing.assert_allclose(tmp[1][0], 0.00098462898029272805)

    hagg = stats.computeHagg(maplayers2['model']['grid'])
    np.testing.assert_allclose(hagg, 60.8999734766)
    parea = stats.computeParea(maplayers2['model']['grid'],
                               probthresh=0.2)
    np.testing.assert_allclose(parea, 48.908522334)

    stats2 = stats.computeStats(maplayers2['model']['grid'], probthresh=0.2, shakefile=shakefile,
                                shakethreshtype='pga', shakethresh=20.,
                                statprobthresh=0.0)
    np.testing.assert_allclose(stats2['Max'], 0.40945028419807472)
    np.testing.assert_allclose(stats2['Median'], 0.00033004360203373138)
    np.testing.assert_allclose(stats2['Std'], 0.04488517841212223)
    np.testing.assert_allclose(stats2['Hagg_0.20g'], 49.70535263249463)
    np.testing.assert_allclose(stats2['Parea_0.20'], 48.908522333887255)

if __name__ == "__main__":
    test_stats_models()
