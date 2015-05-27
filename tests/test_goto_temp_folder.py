from __future__ import print_function
import unittest
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca
from pytraj.utils.context import goto_temp_folder

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]

        def test():
            with goto_temp_folder():
                fa.save("test.nc")
                return mdio.load("test.nc", fa.top)

        fanew = test()
        aa_eq(fanew.xyz, fa.xyz)

        fanew.center()
        fa.center()
        aa_eq(fanew.xyz, fa.xyz)

if __name__ == "__main__":
    unittest.main()
