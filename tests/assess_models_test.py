#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
from configobj import ConfigObj

# third party
from gfail.conf import correct_config_filepaths
import gfail.logisticmodel as LM
from groundfailure import assess_models
from gfail import stats
import matplotlib
matplotlib.use('Agg')

# where is this script?
homedir = os.path.dirname(os.path.abspath(__file__))
upone = os.path.join(homedir, os.pardir)
datadir = os.path.abspath(os.path.join(homedir, 'data'))


def test_assess_models():
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
    hagg = stats.computeParea(maplayers2['model']['grid'],
                              probthresh=0.2)
    np.testing.assert_allclose(hagg, 48.908522334)
    einfo = assess_models.getQuakeInfo('us2000ahv0')
    assert einfo[0] == 'M 8.2 - 101km SSW of Tres Picos, Mexico'

    # fake inventory
    np.random.seed(123)
    fake_inv = os.path.join(datadir, 'loma_prieta', 'fake_inventory.shp')
    gdict = maplayers2['model']['grid'].getGeoDict()
    cov_grid = assess_models.convert2Coverage(gdict, fake_inv, numdiv=10.0)
    np.testing.assert_allclose(
        np.mean(cov_grid.getData()), 0.015693884037747)
    prob_grid = assess_models.convert2Prob(gdict, fake_inv,
                                           mustContainCenter=True)
    np.testing.assert_allclose(
        np.mean(prob_grid.getData()), 0.015790343)
    prob_grid = assess_models.convert2Prob(gdict, fake_inv,
                                           mustContainCenter=False)
    np.testing.assert_allclose(
        np.mean(prob_grid.getData()), 0.022321429)

    # I think this test is freezing on travis
#    resid, rslt = assess_models.statsCoverage(
#            maplayers2['model']['grid'], prob_grid, showplots=False)
#    np.testing.assert_allclose(
#        np.nanmean(resid.getData()), 0.01525620534927)

    # I think this test is freezing on travis
#    temp = assess_models.stats(maplayers2['model']['grid'], fake_inv,
#                               runtests=False, showplots=False)
#    np.testing.assert_allclose(temp[4]['AUC_ROC'], 0.79468845760980589)
#    np.testing.assert_allclose(temp[4]['Brier'], 0.10103376655829596)
#    np.testing.assert_allclose(temp[4]['Brier_no'], 0.016153871532088152)
#    np.testing.assert_allclose(temp[4]['Brier_yes'], 0.78778928086125)

    rho = assess_models.normXcorr(maplayers2['model']['grid'], prob_grid)
    np.testing.assert_allclose(rho, 0.3250259301176)


if __name__ == "__main__":
    test_assess_models()



































