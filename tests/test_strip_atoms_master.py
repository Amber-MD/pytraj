#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestStripAtomsMaster(unittest.TestCase):
    def test_strip_traj_top_fiterator_(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")

        strip_atoms = pt.strip_atoms
        print(strip_atoms(traj, '@CA'))
        print(strip_atoms(traj(), '@CA'))
        print(strip_atoms(traj[:], '@CA'))
        print(strip_atoms(traj.top, '@CA'))


if __name__ == "__main__":
    unittest.main()
