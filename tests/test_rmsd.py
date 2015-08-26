import unittest; import pytraj as pt
import pytraj as pt
from pytraj.base import *
from pytraj import io as mdio
import numpy as np
from pytraj.testing import test_if_having, no_test
from pytraj.utils import assert_almost_equal

TRAJ = TrajectoryIterator(
    filename="./data/md1_prod.Tc5b.x",
    top="./data/Tc5b.top")
cpptraj_rmsd = np.loadtxt(
    "./data/rmsd_to_firstFrame_CA_allres.Tc5b.dat",
    skiprows=1).transpose()[1]


class Test(unittest.TestCase):
    def test_0(self):
        farray = Trajectory()
        farray.top = TRAJ.top
        #print("test_info")
        i = 0
        for frame in TRAJ:
            i += 1
            frame.strip_atoms(mask="!@CA", top=TRAJ.top.copy())
            farray.append(frame.copy())
        assert i == TRAJ.size == TRAJ.n_frames
        assert frame.size == TRAJ.top.n_residues * 3
        farray.top.strip_atoms("!@CA")
        #print("farray.top.n_atoms= ", farray.top.n_atoms)
        assert farray.top.n_atoms == TRAJ.top.n_residues
        farray.top.summary()
        assert farray.size == TRAJ.n_frames
        #print("rmsd to first = ", farray[0].rmsd(farray[1]))
        arr = np.zeros(farray.size)
        #print(cpptraj_rmsd[:10])

        # caculate rmsd to 1st frame
        for i in range(farray.size):
            arr[i] = farray[0].rmsd(farray[i])
        #print(arr[:10])
        np.testing.assert_almost_equal(arr, cpptraj_rmsd, decimal=3)
        #print("Kool, reproduce cpptraj output")

    def test_rmsd_with_mask(self):
        TRAJ = TrajectoryIterator(
            filename="./data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top")
        cpptraj_rmsd = np.loadtxt(
            "./data/rmsd_to_firstFrame_CA_allres.Tc5b.dat",
            skiprows=1).transpose()[1]
        f0 = TRAJ[0]
        arr0 = np.zeros(TRAJ.size)
        arr1 = np.zeros(TRAJ.size)
        mask = "@CA"
        atm = AtomMask(mask)
        TRAJ.top.set_integer_mask(atm)

        for i, frame in enumerate(TRAJ):
            arr0[i] = frame.rmsd(f0, mask=mask, top=TRAJ.top)
            arr1[i] = frame.rmsd(f0, atommask=atm)
        #print(arr0[:10])
        #print(arr1[:10])

        arr2 = pt.rmsd(TRAJ, mask=mask, ref=f0)
        arr3 = pt.rmsd(TRAJ, mask=mask, ref=0)
        np.testing.assert_almost_equal(arr0, cpptraj_rmsd, decimal=3)
        np.testing.assert_almost_equal(arr1, cpptraj_rmsd, decimal=3)
        np.testing.assert_almost_equal(arr2, cpptraj_rmsd, decimal=3)
        np.testing.assert_almost_equal(arr3, cpptraj_rmsd, decimal=3)

    @test_if_having("mdtraj")
    def test_action_rmsd(self):
        # use `mdtraj` for rerefence values
        import mdtraj as md
        traj = TrajectoryIterator(
            filename="./data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top")
        import pytraj.common_actions as pyca
        m_top = md.load_prmtop("./data/Tc5b.top")
        m_traj = md.load_mdcrd("./data/md1_prod.Tc5b.x", m_top)
        m_traj.xyz = m_traj.xyz * 10  # convert `nm` to `Angstrom` unit

        #print("rmsd to first, all atoms")
        arr0 = pt.rmsd(traj, 0)
        arr1 = pt.rmsd(traj, ref=0)
        arr2 = pt.rmsd(traj, )
        a_md0 = md.rmsd(m_traj, m_traj, 0)
        assert_almost_equal(arr0, arr1)
        assert_almost_equal(arr0, arr2)
        assert_almost_equal(arr0, a_md0)

        #print("rmsd to last frame, all atoms")
        arr0 = pt.rmsd(traj, ref=-1)
        arr1 = pt.rmsd(traj, ref=-1)
        a_md = md.rmsd(m_traj, m_traj, -1)
        assert_almost_equal(arr0, arr1)
        assert_almost_equal(arr0, a_md)

        #print("rmsd with mask and indices")
        mask = ":3-18@CA,C"
        atm = traj.top(mask)
        arr0 = pt.rmsd(traj, ref=-1, mask=mask)
        arr1 = pt.rmsd(traj, mask=atm.indices, ref=-1)
        arr2 = pt.rmsd(traj, mask=list(atm.indices), ref=-1)
        arr3 = pt.rmsd(traj, mask=tuple(atm.indices), ref=-1)
        a_md = md.rmsd(m_traj, m_traj, -1, atm.indices)
        #print('arr0', arr0, 'a_md', a_md)
        assert_almost_equal(arr0, a_md)
        assert_almost_equal(arr1, a_md)
        assert_almost_equal(arr2, a_md)
        assert_almost_equal(arr3, a_md)

        from pytraj import Trajectory
        fa = Trajectory(traj)
        arr0 = pt.rmsd(fa, ref=-1, mask=mask)
        arr1 = pt.rmsd(fa, mask=atm.indices, ref=-1)
        arr2 = pt.rmsd(fa, mask=list(atm.indices), ref=-1)
        arr3 = pt.rmsd(fa, mask=tuple(atm.indices), ref=-1)
        a_md = md.rmsd(m_traj, m_traj, -1, atm.indices)
        assert_almost_equal(arr0, a_md)
        assert_almost_equal(arr1, a_md)
        assert_almost_equal(arr2, a_md)
        assert_almost_equal(arr3, a_md)

        # mode = 'cpptraj'
        from pytraj import Trajectory
        fa = Trajectory(traj)
        mask = "!@H="
        atm = fa.top(mask)
        arr0 = pt.rmsd(fa, ref=4, mask=mask)
        a_md = md.rmsd(m_traj, m_traj, 4, atm.indices)
        #print('mode=cpptraj', arr0[1:])
        #print(a_md)

        # exclude 0-th frame for ref
        assert_almost_equal(arr0, a_md)

    def test_list_of_masks(self):
        aa_eq = assert_almost_equal
        traj = TRAJ.copy()
        mask = ['@CA', '@CB', ':3-18@CA,C']
        arr = pt.rmsd(traj, mask=mask)
        for idx, m in enumerate(mask):
            aa_eq(arr[idx], pt.rmsd(traj, mask=m))
            aa_eq(arr[idx], pt.rmsd(traj, mask=traj.top.select(m)))

        mask = ['@CA', '@CB', ':3-18@CA,C', [0, 3, 5]]
        self.assertRaises(ValueError, lambda: pt.rmsd(traj, mask=mask))


if __name__ == "__main__":
    unittest.main()
