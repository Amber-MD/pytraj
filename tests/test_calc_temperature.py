import unittest
from pytraj import *
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        # is_ensemble = False
        topn = "./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7"
        traj = io.iterload("./data/Test_RemdTraj/rem.nc.000", topn)

        # based on existing T
        print(traj.temperatures)


if __name__ == "__main__":
    unittest.main()
