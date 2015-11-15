#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestCreateCRD(unittest.TestCase):

    def test_autoimage_rms_strip(self):
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")

        state = pt.load_batch(traj, '''
        autoimage
        rms
        strip !@CA
        createcrd''')
        state.run()

        # get last datset (DatasetCoordsCRD)
        crd = state.data[-1]
        t0 = traj[:]
        new_traj = t0.autoimage().superpose()['@CA']
        aa_eq(new_traj.xyz, crd.xyz)


if __name__ == "__main__":
    unittest.main()
