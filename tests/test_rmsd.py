import unittest
from pytraj.base import *
from pytraj import io as mdio
import numpy as np
from pytraj.testing import test_if_having
from pytraj.utils import assert_almost_equal

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
        mask = "@CA"
        atm = AtomMask(mask)
        TRAJ.top.set_integer_mask(atm)

        for i, frame in enumerate(TRAJ):
            arr0[i] = frame.rmsd(f0, mask=mask, top=TRAJ.top)
            arr1[i] = frame.rmsd(f0, atommask=atm)
        print(arr0[:10])
        print(arr1[:10])

        arr2 = TRAJ.calc_rmsd(mask, f0) 
        arr3 = TRAJ.calc_rmsd(mask, 0) 
        np.testing.assert_almost_equal(arr0, cpptraj_rmsd, decimal=3)
        np.testing.assert_almost_equal(arr1, cpptraj_rmsd, decimal=3)
        np.testing.assert_almost_equal(arr2, cpptraj_rmsd, decimal=3)
        np.testing.assert_almost_equal(arr3, cpptraj_rmsd, decimal=3)

    @test_if_having("mdtraj")
    def test_action_rmsd(self):
        # use `mdtraj` for rerefence values
        import mdtraj as md
        traj = TrajReadOnly(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
        import pytraj.common_actions as pyca
        m_top = md.load_prmtop("./data/Tc5b.top")
        m_traj = md.load_mdcrd("./data/md1_prod.Tc5b.x", m_top)
        m_traj.xyz = m_traj.xyz * 10 # convert `nm` to `Angstrom` unit

        # rmsd to first, all atoms
        arr0 = traj.calc_rmsd("", 0)
        arr1 = traj.calc_rmsd("", 'first')
        arr2 = traj.calc_rmsd()
        a_md0 = md.rmsd(m_traj, m_traj, 0)
        assert_almost_equal(arr0, arr1)
        assert_almost_equal(arr0, arr2)
        assert_almost_equal(arr0, a_md0)

        # rmsd to last frame, all atoms
        arr0 = traj.calc_rmsd(ref='last')
        arr1 = traj.calc_rmsd(ref=-1)
        a_md = md.rmsd(m_traj, m_traj, -1)
        assert_almost_equal(arr0, arr1)
        assert_almost_equal(arr0, a_md)

        # rmsd with mask and indices
        mask = ":3-18@CA, C"
        atm = traj.top(mask)
        arr0 = traj.calc_rmsd(mask, ref='last')
        arr1 = traj.calc_rmsd(atm.indices, ref='last')
        arr2 = traj.calc_rmsd(list(atm.indices), ref='last')
        arr3 = traj.calc_rmsd(tuple(atm.indices), ref='last')
        a_md = md.rmsd(m_traj, m_traj, -1, atm.indices)
        assert_almost_equal(arr0, a_md)
        assert_almost_equal(arr1, a_md)
        assert_almost_equal(arr2, a_md)
        assert_almost_equal(arr3, a_md)

if __name__ == "__main__":
    unittest.main()
