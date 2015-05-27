from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from pytraj import AtomMask

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0 = traj[0]
        f0s = f0.copy()
        xyzf4 = f0.xyz.astype('f4').copy()

        # cast from np f4
        f0[:] = xyzf4 # cast from f4 to f8
        aa_eq(f0.xyz, f0s.xyz)

        # an integer
        f0[0, 0] = 1
        assert f0[0, 0] == 1.

        # dummy cast
        f0[0] = 1
        aa_eq(f0[0], [1., 1., 1.])

        indices = traj.top("@CA").indices
        newxyz = xyzf4[indices]
        newxyz += 1.

        f0[indices] = newxyz
        aa_eq(f0[indices], newxyz)

        atm = AtomMask(indices)
        aa_eq(f0[atm], f0[indices])

        newxyz += 2.
        f0[atm] = newxyz
        aa_eq(f0[atm], f0[indices])
        aa_eq(f0[indices], newxyz)

        f0[atm] = newxyz.tolist()
        aa_eq(f0[atm], f0[indices])
        aa_eq(f0[indices], newxyz)

        f0[traj.top('@CA')] = newxyz.tolist()
        aa_eq(f0[traj.top('@CA')], newxyz)


if __name__ == "__main__":
    unittest.main()
