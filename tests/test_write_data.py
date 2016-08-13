#!/usr/bin/env python

from __future__ import print_function
import os
import unittest
import numpy as np
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestWriteData(unittest.TestCase):

    def test_write_data(self):
        traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")
        rmsd_data = pt.rmsd(traj)

        try:
            os.remove("test.agr")
            os.remove("test.gnu")
            os.remove("test.dx")
        except OSError:
            pass

        pt.io.write_data("test.agr", rmsd_data)
        pt.io.write_data("test.gnu", rmsd_data)

        assert os.path.exists("test.agr")
        assert os.path.exists("test.gnu")

        # not supported yet
        # data_3d = pt.gist(traj, options='name gist')['gist[gO]']
        # pt.io.write_data("test.dx", data_3d)
        # assert os.path.exists("test.dx")

if __name__ == "__main__":
    unittest.main()
