#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
from configobj import ConfigObj

# third party
from mapio.gmt import GMTGrid  # REPLACE THIS

from gfail.conf import correct_config_filepaths
from gfail import (
    Zhu2015Model,
    Nowicki2014Model,
    Zhu2017Model,
    Zhu2017ModelCoastal,
    Jessee2018Model,
)
from gfail import godt2008


homedir = os.path.dirname(os.path.abspath(__file__))  # find where this script is
upone = os.path.join(homedir, os.pardir)
datadir = os.path.abspath(os.path.join(homedir, "data"))

changetarget = False  # Turn to True if need to recompute target data


def test_zhu2015():
    conf_file = os.path.join(upone, "defaultconfigfiles", "models", "zhu_2015.ini")
    conf = ConfigObj(conf_file)
    data_path = os.path.join(datadir, "loma_prieta", "model_inputs")
    # Check slopefile trimming
    conf["zhu_2015"]["slopefile"] = "global_gted_maxslope_30c.flt"
    conf = correct_config_filepaths(data_path, conf)
    # Run with divfactor of 1
    conf["zhu_2015"]["divfactor"] = "1."

    shakefile = os.path.join(datadir, "loma_prieta", "grid.xml")
    lm = Zhu2015Model(
        shakefile,
        conf["zhu_2015"],
        uncertfile=None,
        bounds=None,
        trimfile=None,
        saveinputs=True,
    )
    maplayers = lm.calculate()

    pgrid = maplayers["model"]["grid"]
    test_data = pgrid.getData()

    if changetarget:
        # To change target data:
        pgrd = GMTGrid(pgrid.getData(), pgrid.getGeoDict())
        pgrd.save(os.path.join(datadir, "loma_prieta", "targets", "zhu2015.grd"))

    # Load target
    target_file = os.path.join(datadir, "loma_prieta", "targets", "zhu2015.grd")
    target_grid = GMTGrid.load(target_file)
    target_data = target_grid.getData()
    # Assert
    np.testing.assert_allclose(target_data, test_data, rtol=1e-3)


def test_zhu_2017_general():
    conf_file = os.path.join(
        upone, "defaultconfigfiles", "models", "zhu_2017_general.ini"
    )
    conf = ConfigObj(conf_file)
    data_path = os.path.join(datadir, "loma_prieta", "model_inputs")
    conf = correct_config_filepaths(data_path, conf)
    # Run with divfactor of 4
    conf["zhu_2017_general"]["divfactor"] = "4."
    shakefile = os.path.join(datadir, "loma_prieta", "grid.xml")
    uncertainty_file = os.path.join(datadir, "loma_prieta", "uncertainty.xml")
    lm = Zhu2017Model(
        shakefile,
        conf["zhu_2017_general"],
        uncertfile=uncertainty_file,
        bounds=None,
        trimfile=None,
        saveinputs=True,
    )
    maplayers = lm.calculate()

    pgrid = maplayers["model"]["grid"]
    test_data = pgrid.getData()

    if changetarget:
        # To change target data:
        pgrd = GMTGrid(pgrid.getData(), pgrid.getGeoDict())
        pgrd.save(
            os.path.join(datadir, "loma_prieta", "targets", "zhu2017_general_div4.grd")
        )

    # Load target
    target_file = os.path.join(
        datadir, "loma_prieta", "targets", "zhu2017_general_div4.grd"
    )
    target_grid = GMTGrid.load(target_file)
    target_data = target_grid.getData()

    # Assert
    np.testing.assert_allclose(target_data, test_data, rtol=1e-3)


def test_zhu_2017_coastal():
    conf_file = os.path.join(
        upone, "defaultconfigfiles", "models", "zhu_2017_coastal.ini"
    )
    conf = ConfigObj(conf_file)
    data_path = os.path.join(datadir, "loma_prieta", "model_inputs")
    conf = correct_config_filepaths(data_path, conf)
    # Run with divfactor of 1
    conf["zhu_2017_coastal"]["divfactor"] = "1."
    shakefile = os.path.join(datadir, "loma_prieta", "grid.xml")
    lm = Zhu2017ModelCoastal(
        shakefile,
        conf["zhu_2017_coastal"],
        uncertfile=None,
        bounds=None,
        trimfile=None,
        saveinputs=True,
    )
    maplayers = lm.calculate()

    pgrid = maplayers["model"]["grid"]
    test_data = pgrid.getData()

    if changetarget:
        # To change target data:
        pgrd = GMTGrid(pgrid.getData(), pgrid.getGeoDict())
        pgrd.save(
            os.path.join(datadir, "loma_prieta", "targets", "zhu2017_coastal.grd")
        )

    # Load target
    target_file = os.path.join(datadir, "loma_prieta", "targets", "zhu2017_coastal.grd")
    target_grid = GMTGrid.load(target_file)
    target_data = target_grid.getData()

    # Assert
    np.testing.assert_allclose(target_data, test_data, rtol=1e-3)


def test_nowicki_2014_global():
    conf_file = os.path.join(
        upone, "defaultconfigfiles", "models", "nowicki_2014_global.ini"
    )
    conf = ConfigObj(conf_file)
    data_path = os.path.join(datadir, "loma_prieta", "model_inputs")
    conf = correct_config_filepaths(data_path, conf)
    # Run with divfactor of 1
    conf["nowicki_2014_global"]["divfactor"] = "1."
    shakefile = os.path.join(datadir, "loma_prieta", "grid.xml")
    lm = Nowicki2014Model(
        shakefile,
        conf["nowicki_2014_global"],
        uncertfile=None,
        bounds=None,
        trimfile=None,
        saveinputs=True,
    )
    maplayers = lm.calculate()

    pgrid = maplayers["model"]["grid"]
    test_data = pgrid.getData()

    if changetarget:
        # To change target data:
        pgrd = GMTGrid(pgrid.getData(), pgrid.getGeoDict())
        pgrd.save(
            os.path.join(datadir, "loma_prieta", "targets", "nowicki_2014_global.grd")
        )

    # Load target
    target_file = os.path.join(
        datadir, "loma_prieta", "targets", "nowicki_2014_global.grd"
    )
    target_grid = GMTGrid.load(target_file)
    target_data = target_grid.getData()

    # Assert
    np.testing.assert_allclose(target_data, test_data, rtol=1e-3)


def test_jessee_2018():
    conf_file = os.path.join(upone, "defaultconfigfiles", "models", "jessee_2018.ini")
    conf = ConfigObj(conf_file)
    data_path = os.path.join(datadir, "loma_prieta", "model_inputs")
    conf = correct_config_filepaths(data_path, conf)
    shakefile = os.path.join(datadir, "loma_prieta", "grid.xml")
    uncertainty_file = os.path.join(datadir, "loma_prieta", "uncertainty.xml")
    lm = Jessee2018Model(
        shakefile,
        conf["jessee_2018"],
        uncertfile=uncertainty_file,
        bounds=None,
        trimfile=None,
        saveinputs=True,
    )
    maplayers = lm.calculate()

    pgrid = maplayers["model"]["grid"]
    stdgrid = maplayers["std"]["grid"]
    test_data = pgrid.getData()
    test_data_std = stdgrid.getData()

    if changetarget:
        # To change target data:
        pgrd = GMTGrid(pgrid.getData(), pgrid.getGeoDict())
        pgrd.save(os.path.join(datadir, "loma_prieta", "targets", "jessee_2018.grd"))
        stdgrd = GMTGrid(stdgrid.getData(), stdgrid.getGeoDict())
        stdgrd.save(
            os.path.join(datadir, "loma_prieta", "targets", "jessee_2018_std.grd")
        )

    # Load target
    target_file = os.path.join(datadir, "loma_prieta", "targets", "jessee_2018.grd")
    target_grid = GMTGrid.load(target_file)
    target_data = target_grid.getData()
    std_file = os.path.join(datadir, "loma_prieta", "targets", "jessee_2018_std.grd")
    target_grid_std = GMTGrid.load(std_file)
    target_data_std = target_grid_std.getData()

    # Assert
    np.testing.assert_allclose(target_data, test_data, rtol=1e-3)
    np.testing.assert_allclose(target_data_std, test_data_std, rtol=1e-3)


def test_godt_2008():
    conf_file = os.path.join(upone, "defaultconfigfiles", "models", "godt_2008.ini")
    conf = ConfigObj(conf_file)
    conf["godt_2008"]["divfactor"] = "1."
    data_path = os.path.join(datadir, "loma_prieta", "model_inputs")
    conf = correct_config_filepaths(data_path, conf)
    shakefile = os.path.join(datadir, "loma_prieta", "grid.xml")
    maplayers = godt2008(shakefile, conf)

    pgrid = maplayers["model"]["grid"]
    test_data = pgrid.getData()

    if changetarget:
        # To change target data:
        pgrd = GMTGrid(pgrid.getData(), pgrid.getGeoDict())
        pgrd.save(os.path.join(datadir, "loma_prieta", "targets", "godt_2008.grd"))

    # Load target
    target_file = os.path.join(datadir, "loma_prieta", "targets", "godt_2008.grd")
    target_grid = GMTGrid.load(target_file)
    target_data = target_grid.getData()

    # Assert
    np.testing.assert_allclose(target_data, test_data, rtol=1e-3)


if __name__ == "__main__":
    test_zhu2015()
    test_zhu_2017_general()
    test_zhu_2017_coastal()
    test_nowicki_2014_global()
    test_jessee_2018()
    test_godt_2008()
