from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    def test_0(self):
        from pytraj import _auto_funcs as af
        from pytraj import info

        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        # wrong result with
        # d0 = af.calc_molsurf(traj).values
        # d0 = af.calc_molsurf(traj, dtype='ndarray')
        # don't need to fix now since `_auto_funcs` is private method,
        # just for testing purpose

        d0 = af.calc_molsurf(traj)
        d1 = traj.calc_molsurf(dtype='dataset').to_ndarray()
        print(d0.values)
        print(d1)
        aa_eq(d0.values, d1)

if __name__ == "__main__":
    unittest.main()
