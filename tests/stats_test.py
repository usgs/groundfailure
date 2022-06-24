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

matplotlib.use("Agg")

# where is this script?
homedir = os.path.dirname(os.path.abspath(__file__))
upone = os.path.join(homedir, os.pardir)
datadir = os.path.abspath(os.path.join(homedir, "data"))


def test_stats_models():
    conf_file = os.path.join(upone, "oldconfigfiles", "models", "zhu_2015.ini")
    conf = ConfigObj(conf_file)
    data_path = os.path.join(datadir, "loma_prieta", "model_inputs")
    # Check slopefile trimming
    conf["zhu_2015"]["slopefile"] = "global_gted_maxslope_30c.flt"
    conf = correct_config_filepaths(data_path, conf)
    # Run with divfactor of 1
    conf["zhu_2015"]["divfactor"] = "1."

    shakefile = os.path.join(datadir, "loma_prieta", "grid.xml")
    lm = LM.LogisticModel(shakefile, conf, saveinputs=True)
    # maplayers1 = lm.calculate()
    conf_file = os.path.join(upone, "oldconfigfiles", "models", "zhu_2017_coastal.ini")
    conf = ConfigObj(conf_file)
    data_path = os.path.join(datadir, "loma_prieta", "model_inputs")
    conf["zhu_2017_coastal"]["slopefile"] = "global_gted_maxslope_30c.flt"
    conf = correct_config_filepaths(data_path, conf)
    # Run with divfactor of 1
    conf["zhu_2017_coastal"]["divfactor"] = "1."
    shakefile = os.path.join(datadir, "loma_prieta", "grid.xml")
    lm = LM.LogisticModel(shakefile, conf, saveinputs=True)
    maplayers2 = lm.calculate()
    # Change shakemap name so that it doesn't stomp on the other
    maplayers2["model"]["description"]["shakemap"] = "19891018000415_ver2"

    hagg = stats.computeHagg(maplayers2["model"]["grid"])
    np.testing.assert_allclose(hagg["hagg_0.00g"], 5.155065, atol=0.001)

    stats2 = stats.computeStats(
        maplayers2["model"]["grid"],
        shakefile=shakefile,
        shakethreshtype="pga",
        shakethresh=20.0,
        probthresh=0.0,
    )
    np.testing.assert_allclose(stats2["Max"], 0.4105819792026343, atol=0.001)
    np.testing.assert_allclose(stats2["Median"], 0.34855563636356035, rtol=0.01)
    np.testing.assert_allclose(stats2["Std"], 0.04195, atol=0.001)
    np.testing.assert_allclose(stats2["hagg_0.20g"], 2.854704, atol=0.001)


if __name__ == "__main__":
    test_stats_models()
