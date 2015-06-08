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
        traj = mdio.load_sample_data("tz2")[:]
        traj.autoimage()
        traj.rmsfit(mask=':1-13')
        d = pyca.calc_grid(traj, " 20 0.5 20 0.5 20 0.5 :WAT@O")
        d[0].save("./output/test_grid.dx")
        print (d[0])

        # iterator
        d = pyca.calc_grid(traj(), " 20 0.5 20 0.5 20 0.5 :WAT@O", top=traj.top)
        d[0].save("./output/test_grid2.dx")
        print (d[0])

if __name__ == "__main__":
    unittest.main()
