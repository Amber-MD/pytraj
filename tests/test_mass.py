import unittest
from array import array
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        top = traj.top

        mlist = []
        for atom in top:
            mlist.append(atom.mass)
        mlist = array('d', mlist)

        assert_almost_equal(top.mass, mlist)


if __name__ == "__main__":
    unittest.main()
