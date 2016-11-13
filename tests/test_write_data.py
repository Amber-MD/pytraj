#!/usr/bin/env python

from __future__ import print_function
import os
import unittest
import pytraj as pt
from pytraj.testing import cpptraj_test_dir, tempfolder


def test_write_data():
    traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")
    rmsd_data = pt.rmsd(traj)

    with tempfolder():
        try:
            os.remove("test.agr")
            os.remove("test.gnu")
            os.remove("test.dx")
            os.remove("test.ccp4")
        except OSError:
            pass

        pt.io.write_data("test.agr", rmsd_data)
        pt.io.write_data("test.gnu", rmsd_data)

        assert os.path.exists("test.agr")
        assert os.path.exists("test.gnu")


        saved_dx_file = cpptraj_test_dir + '/Test_CCP4/fav8.dx.save'
        output_file = 'test.ccp4'
        pt.datafiles.convert(saved_dx_file, output_file)
        assert os.path.exists("test.ccp4")
