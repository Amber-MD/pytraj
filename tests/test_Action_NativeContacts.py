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
        # TODO : assert
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        ref = mdio.load("data/Tc5b.nat.crd", traj.top)
        dslist = pyca.native_contacts(traj, top=traj.top, ref=ref)
        #print(dslist)

        dslist = pyca.native_contacts(traj, top=traj.top)
        #print(dslist)


if __name__ == "__main__":
    unittest.main()
