import unittest
from array import array
import numpy as np
import pytraj as pt
from pytraj.testing import aa_eq
from pytraj import Frame

SMALL = 1E-6

N_ATOMS = 10
FRAME = Frame(N_ATOMS)
arr = np.arange(3 * N_ATOMS)
FRAME.xyz[:] = arr.reshape(N_ATOMS, 3)
FRAME_orig = FRAME.copy()


class TestFrame(unittest.TestCase):

    def test_fit(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        trajnew = pt.iterload("./data/md1_prod.fit_to_first.Tc5b.x",
                              "./data/Tc5b.top")

        # make sure 0-th frame does not change
        frame0 = traj[0]
        frame0new = trajnew[0]
        aa_eq(frame0.xyz, frame0new.xyz)

        frame1 = traj[1]
        frame1new = trajnew[1]

        # try do-fitting from Python
        # not right yet
        rmsd, mat, v1, v2 = frame1.rmsd(frame0, get_mvv=True)
        frame1._trans_rot_trans(v1, mat, v2)
        assert frame1.rmsd(frame1new) < 1E-3

    def test_1(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        trajnew = pt.iterload("./data/md1_prod.fit_to_first.Tc5b.x",
                              "./data/Tc5b.top")
        frame0 = traj[0]
        assert frame0[0, 2] == frame0.atom(0)[2]
        aa_eq(frame0[0, :], frame0.atom(0))


if __name__ == "__main__":
    unittest.main()
