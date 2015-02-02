import unittest
from pytraj.base import *
from pytraj import io as mdio
import numpy as np

TRAJ = TrajReadOnly(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
cpptraj_rmsd = np.loadtxt("./data/rmsd_to_firstFrame_CA_allres.Tc5b.dat", skiprows=1).transpose()[1]

class Test(unittest.TestCase):
    def test_0(self):
        farray = FrameArray()
        farray.top = TRAJ.top
        print("test_info")
        i = 0
        for frame in TRAJ:
            i +=1
            frame.strip_atoms(mask="!@CA", top=TRAJ.top.copy())
            farray.append(frame)
        assert i == TRAJ.size == TRAJ.max_frames
        assert frame.size == TRAJ.top.n_res * 3
        farray.top.strip_atoms("!@CA")
        print("farray.top.n_atoms= ", farray.top.n_atoms)
        assert farray.top.n_atoms == TRAJ.top.n_res 
        farray.top.summary()
        assert farray.size == TRAJ.max_frames
        print("rmsd to first = ", farray[0].rmsd(farray[1]))
        arr = np.zeros(farray.size)
        print(cpptraj_rmsd[:10])

        # caculate rmsd to 1st frame
        for i in range(farray.size):
            arr[i] = farray[0].rmsd(farray[i])
        print(arr[:10])
        np.testing.assert_almost_equal(arr, cpptraj_rmsd, decimal=3)
        print("Kool, reproduce cpptraj output")

    def test_rmsd_with_mask(self):
        TRAJ = TrajReadOnly(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
        cpptraj_rmsd = np.loadtxt("./data/rmsd_to_firstFrame_CA_allres.Tc5b.dat", skiprows=1).transpose()[1]
        f0 = TRAJ[0]
        arr0 = np.zeros(TRAJ.size)
        arr1 = np.zeros(TRAJ.size)
        atm = AtomMask("@CA")
        TRAJ.top.set_integer_mask(atm)

        for i, frame in enumerate(TRAJ):
            arr0[i] = frame.rmsd(f0, mask="@CA", top=TRAJ.top)
            arr1[i] = frame.rmsd(f0, atommask=atm)
        print(arr0[:10])
        print(arr1[:10])
        np.testing.assert_almost_equal(arr0, cpptraj_rmsd, decimal=3)
        np.testing.assert_almost_equal(arr1, cpptraj_rmsd, decimal=3)

if __name__ == "__main__":
    unittest.main()
