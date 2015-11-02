#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq

class TestSuperDispatch(unittest.TestCase):
    def test_super_dispatch(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")

        funclist = [pt.radgyr, pt.molsurf]
        for func in funclist:
            mask = '@CA'
            atom_indices = pt.select_atoms(traj.top, mask)
            # mask
            aa_eq(func(traj, mask=mask),
                  func(traj, mask=atom_indices))

            # frame_indices with mask
            frame_indices = [0, 5, 8]
            aa_eq(func(traj[frame_indices], mask=mask),
                  func(traj, mask=atom_indices, frame_indices=frame_indices))


if __name__ == "__main__":
    unittest.main()
