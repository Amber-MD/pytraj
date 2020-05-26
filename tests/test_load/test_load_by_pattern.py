from __future__ import print_function
import pytraj as pt
from utils import fn
import unittest
from pytraj.testing import aa_eq
import pytest


class Test(unittest.TestCase):
    def test_load_from_list(self):
        from glob import glob
        pattern = fn('Test_RemdTraj/rem.nc.*')
        flist = sorted(glob(pattern))
        top = glob(fn('Test_RemdTraj/ala*parm7'))[0]
        traj0 = pt.iterload(flist, top)
        traj1 = pt.iterload(pattern, top)
        aa_eq(traj0.xyz, traj1.xyz)

        # raise if not find files
        with pytest.raises(ValueError):
            pt.iterload("./data/xyz_cool*.x", traj0.top)


if __name__ == "__main__":
    unittest.main()
