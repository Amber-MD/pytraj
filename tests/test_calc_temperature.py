import unittest
from pytraj import *


class Test(unittest.TestCase):

    def test_0(self):
        topn = "./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7"
        traj = io.iterload("./data/Test_RemdTraj/rem.nc.000", topn)

        # based on existing T


if __name__ == "__main__":
    unittest.main()
