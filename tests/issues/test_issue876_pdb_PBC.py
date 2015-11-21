#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestPBCFromPDB(unittest.TestCase):

    def test_pbc(self):
        traj = pt.load('data/small_pbc.pdb')
        aa_eq(traj.unitcells[0], [51.263, 51.263, 51.263, 90.00, 90.00, 90.00])
        assert traj.top.has_box() == True, 'Topology must has box'
        assert traj.top.box.type == 'ortho', 'must be ortho box'


if __name__ == "__main__":
    unittest.main()
