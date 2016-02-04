import unittest
import pytraj as pt
from pytraj.base import *
from pytraj.testing import aa_eq
from pytraj import io as mdio
import numpy as np
from array import array


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        assert frame0._buffer2d[0, 0] == frame0[0, 0]
        assert frame0._buffer2d[0, 1] == frame0[0, 1]
        assert frame0._buffer2d[0, 2] == frame0[0, 2]
        assert frame0._buffer2d[1, 0] == frame0[1, 0]
        assert frame0._buffer2d[19, 0] == frame0[19, 0]
        frame0._buffer2d[19, 0] = 1000.

        # make sure changing _buffer2d will update frame0.xyz too
        assert frame0._buffer2d[19, 0] == frame0[19, 0]
        arr0 = np.asarray(frame0._buffer2d)
        arr0[19] = [200, 300, 400.]
        assert frame0._buffer2d[19, 0] == frame0[19, 0] == arr0[19, 0]

        frame1 = frame0.copy()
        frame1.strip(traj.top("!@CA"))

        _f = pt.Frame(frame0, traj.top("@CA"))
        CA_2 = _f[2, :]
        aa_eq(CA_2, frame1[2, :])
        assert frame0._buffer1d.shape == (912, )
        assert frame0._buffer2d.shape == (304, 3)
        assert frame0._buffer1d.is_c_contig() == True
        frame0._buffer2d[1:3, 0] = array('d', [1., 2.])
        aa_eq(frame0[1:3, 0], array('d', frame0._buffer2d[1:3, 0]))
        aa_eq(frame0._buffer2d[1:3, 0], array('d', frame0._buffer2d[1:3, 0]))

    def test_1(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        frame0[1:5, 2] = list(range(100, 104))
        frame0[0, :] = list(range(3))

    def test_magic_methods(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        frame1 = frame0.copy()
        frame1 += frame1
        aa_eq(2 * frame0.xyz[0], frame1.xyz[0])


if __name__ == "__main__":
    unittest.main()
