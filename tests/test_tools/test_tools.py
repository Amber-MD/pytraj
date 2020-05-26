#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
import pytest


class TestTools(unittest.TestCase):
    def test_tools(self):
        traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))
        aa_eq(pt.tools.as_2darray(traj), pt.tools.as_2darray(traj.xyz))

        # as_2darray
        assert pt.tools.as_2darray(traj).ndim == 2, 'ndim must be 2'

        #
        with pytest.raises(ValueError):
            pt.tools.rmsd_1darray([3, 2], [[2, 3], [4, 6]])
        with pytest.raises(ValueError):
            pt.tools.rmsd_1darray([[2, 3], [4, 6]], [2, 3])

        # rmsd
        with pytest.raises(ValueError):
            pt.tools.rmsd([[2, 3], [4, 6]], [2, 3])
        with pytest.raises(ValueError):
            pt.tools.rmsd([[2, 3]], [2, 3], flatten=False)

        #
        for frame in pt.tools.split_traj_by_residues(traj, 0, 12):
            pass

        assert not [
            name for name in pt.tools.dir_(frame) if name.startswith('_')
        ]


if __name__ == "__main__":
    unittest.main()
