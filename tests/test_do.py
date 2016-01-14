#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestDo(unittest.TestCase):

    def test_pytraj_do(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")

        ref0 = pt.autoimage(traj[0], top=traj.top)
        ref1 = pt.autoimage(traj[1], top=traj.top)

        data = pt.tools.dict_to_ndarray(pt.compute(
            ['autoimage', 'radgyr nomax @CA', 'rms refindex 0',
             'rms refindex 1'],
            traj,
            ref=[ref0, ref1]))

        t0 = pt.autoimage(traj[:])
        aa_eq(pt.radgyr(t0, '@CA'), data[0])
        aa_eq(pt.rmsd(t0, ref=ref0), data[1])
        aa_eq(pt.rmsd(t0, ref=ref1), data[2])

        # pytraj's method
        aa_eq(pt.compute(pt.radgyr, t0, '@CA'), data[0])


if __name__ == "__main__":
    unittest.main()
