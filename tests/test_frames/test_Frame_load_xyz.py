from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from pytraj import Frame
from pytraj.testing import aa_eq


class Test(unittest.TestCase):

    def test_xyz(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        frame = Frame()
        frame.append_xyz(traj[0].xyz)
        aa_eq(frame.xyz, traj[0].xyz)
        aa_eq(frame.xyz.flatten(), traj[0].xyz.flatten())
        aa_eq(np.array(frame._buffer1d), traj[0].xyz.flatten())


if __name__ == "__main__":
    unittest.main()
