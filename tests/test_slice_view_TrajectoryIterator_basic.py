from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    def test_0(self):
        # Just want to test if getting segfault
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        #fa = fa.copy()
        fa._fast_strip_atoms('@CA')
        fa.strip_atoms('@CA')
        fa = traj[:]
        #fa._fast_strip_atoms('!@CA')
        fa.strip_atoms('!@CA')
        fa = traj[:]
        #fa._fast_strip_atoms('@CA,CB')
        fa.strip_atoms('@CA,CB')
        fa = traj[:]
        #fa._fast_strip_atoms('@CA,CG')
        fa.strip_atoms('@CA,CG')
        fa = traj[:]

if __name__ == "__main__":
    unittest.main()
