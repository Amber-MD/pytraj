import unittest
from array import array
import numpy as np
#from numpy.testing import assert_almost_equal
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj import Frame
from pytraj.base import *
from pytraj import io as mdio

SMALL = 1E-6

N_ATOMS = 10
FRAME = Frame(N_ATOMS)
arr = np.arange(3 * N_ATOMS)
FRAME.set_from_crd(arr)
FRAME_orig = FRAME.copy()


class TestFrame(unittest.TestCase):
    def test_fit(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        trajnew = mdio.iterload(
            "./data/md1_prod.fit_to_first.Tc5b.x", "./data/Tc5b.top")

        # make sure 0-th frame does not change
        frame0 = traj[0]
        frame0new = trajnew[0]
        assert (frame0.coords == frame0new.coords) == True
        #print(frame0[:10])
        #print(frame0new[:10])

        frame1 = traj[1]
        frame1new = trajnew[1]
        #print(frame1[:10])
        #print(frame1new[:10])
        assert (frame1.coords == frame1new.coords) == False

        # try do-fitting from Python
        # not right yet
        rmsd, mat, v1, v2 = frame1.rmsd(frame0, get_mvv=True)
        #print(rmsd)
        frame1.trans_rot_trans(v1, mat, v2)
        #print(frame1[:10])
        #print(frame1.rmsd_nofit(frame1new))
        assert frame1.rmsd(frame1new) < 1E-3

    def test_1(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        trajnew = mdio.iterload(
            "./data/md1_prod.fit_to_first.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        assert frame0[0, 2] == frame0.atoms(0)[2]
        assert_almost_equal(frame0[0, :], frame0.atoms(0))
        #print(frame0[0, :])
        #print(frame0.atoms(0))

        #print(frame0[:, 2])
        framesub = frame0.get_subframe("@CA", traj.top)
        assert framesub.n_atoms == 20
        #print(framesub[19, :])


if __name__ == "__main__":
    unittest.main()
