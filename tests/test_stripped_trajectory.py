#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestStrippedTrajectory(unittest.TestCase):

    def test_stripped_trajectory(self):
        traj_on_disk = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")

        straj = traj_on_disk.autoimage().superpose('@CA').strip(":WAT")
        traj_on_mem_strip = traj_on_disk['!:WAT']

        aa_eq(traj_on_disk['!:WAT'].xyz, straj.xyz)
        aa_eq(traj_on_mem_strip.xyz, straj.xyz)
        aa_eq(traj_on_mem_strip[0].xyz, straj[0].xyz)
        aa_eq(traj_on_mem_strip[:].xyz, straj[:].xyz)
        aa_eq(traj_on_mem_strip[:3].xyz, straj[:3].xyz)
        aa_eq(traj_on_mem_strip[3:10].xyz, straj[3:10].xyz)

        # save
        fn = 'output/test.nc'
        straj.save(fn, overwrite=True)
        aa_eq(straj.xyz, pt.load(fn, straj.top).xyz)

if __name__ == "__main__":
    unittest.main()
