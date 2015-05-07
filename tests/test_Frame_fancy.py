import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame = traj[0]
        frame.set_top(traj.top)

        # print 2D array for CA atoms
        arr0 = frame['@CA']
        # flat list
        arr1 = [item for sublist in arr0 for item in sublist]

        # make sure we do it right
        topCA = traj.top.copy()
        topCA.strip_atoms("!@CA")
        trajCA = mdio.iterload("./data/stripAllButCA.Tc5b.x", topCA)
        print (arr1)
        print (trajCA[0].coords)
        assert_almost_equal(arr1, trajCA[0].coords)

if __name__ == "__main__":
    unittest.main()
