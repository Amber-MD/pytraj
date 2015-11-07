#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestFrameIndices(unittest.TestCase):
    def setUp(self):
        self.traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
    def test_frame_indices_from_yield(self):
        '''extensive and seperated testsing
        '''
        traj = self.traj

        def gen_int():
            for i in range(0, 10, 2):
                yield i

        for idx, frame in enumerate(pt.iterframe(traj, frame_indices=gen_int())):
            pass

        assert idx == len(range(0, 10, 2)) - 1
        aa_eq(frame.xyz, traj[8].xyz)

    def test_frame_indices_for_function(self):
        traj = self.traj
        funclist = [pt.radgyr, ]
        frame_indices = [0, 5, 2]

        for func in funclist:
            aa_eq(func(traj, frame_indices=frame_indices),
                  func(traj[frame_indices]))



if __name__ == "__main__":
    unittest.main()
