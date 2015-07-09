import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
import numpy as np
from array import array


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        #arr0 = np.asarray(frame0)
        print(frame0[0, :])
        assert frame0.buffer2d[0, 0] == frame0[0, 0]
        assert frame0.buffer2d[0, 1] == frame0[0, 1]
        assert frame0.buffer2d[0, 2] == frame0[0, 2]
        assert frame0.buffer2d[1, 0] == frame0[1, 0]
        assert frame0.buffer2d[19, 0] == frame0[19, 0]
        frame0.buffer2d[19, 0] = 1000.

        # make sure changing buffer2d will update frame0.coords too
        assert frame0.buffer2d[19, 0] == frame0[19, 0]
        arr0 = np.asarray(frame0.buffer2d)
        print(arr0.shape)
        arr0[19] = [200, 300, 400.]
        assert frame0.buffer2d[19, 0] == frame0[19, 0] == arr0[19, 0]

        print("try to strip atoms")
        frame1 = frame0.copy()
        frame1.strip_atoms("!@CA", traj.top)

        frame2 = frame0.get_subframe("@CA", traj.top)
        CA_2 = frame2[2, :]
        print(CA_2)

        # this one will mess up with python2.6
        # TODO: fix python2.6: no way
        # got [0, 0, 0] for CA_3
        # (must be [-16.761   7.688  -3.829])
        CA_3 = frame0.get_subframe("@CA", traj.top)[2, :]
        print(CA_3)

        print(frame1[2, :])
        print("compare get_subframe to strip_atoms")
        assert_almost_equal(frame2[2, :], frame1[2, :])
        assert_almost_equal(frame1.coords[6:9], frame1[2, :])
        assert_almost_equal(frame2.coords, frame1.coords)

if __name__ == "__main__":
    unittest.main()
