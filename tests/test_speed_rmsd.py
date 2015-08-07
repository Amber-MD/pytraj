from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
from pytraj.utils import Timer
import pytraj.common_actions as pyca

tip3p_dir = "./data/nogit/tip3p/"


class Test(unittest.TestCase):
    @test_if_path_exists(tip3p_dir)
    def test_0(self):
        # iterload
        @Timer()
        def test_rmsd_pytraj_mode():
            fname = tip3p_dir + "/md.trj"
            topname = tip3p_dir + "/tc5bwat.top"
            traj = mdio.iterload(fname, topname)
            traj.rmsd(ref=10, mode='pytraj')

        @Timer()
        def test_rmsd_cpptraj_mode():
            fname = tip3p_dir + "/md.trj"
            topname = tip3p_dir + "/tc5bwat.top"
            traj = mdio.iterload(fname, topname)
            traj.rmsd(ref=10, mode='cpptraj')

        from pytraj.utils import _import
        has_mdtraj, md = _import("mdtraj")

        if has_mdtraj:

            @Timer()
            def test_rmsd_mdtraj():
                fname = tip3p_dir + "/md.trj"
                topname = tip3p_dir + "/tc5bwat.top"
                mtop = md.load_prmtop(topname)
                m_traj = md.load_netcdf(fname, top=mtop)
                md.rmsd(m_traj, m_traj, 10)

        print("test_rmsd_pytraj_mode")
        test_rmsd_pytraj_mode()
        print("test_rmsd_cpptraj_mode")
        test_rmsd_cpptraj_mode()
        if has_mdtraj:
            print("test_rmsd_mdtraj")
            test_rmsd_mdtraj()

    @test_if_having("mdtraj")
    def test_1(self):
        import mdtraj as md
        # in memory
        fname = tip3p_dir + "/md.trj"
        topname = tip3p_dir + "/tc5bwat.top"
        traj = mdio.load(fname, topname, indices=range(1000))
        print(traj)

        @Timer()
        def test_rmsd_cpptraj_mode():
            traj.rmsd(ref=10, mode='cpptraj')

        print("test_rmsd_cpptraj_mode")
        test_rmsd_cpptraj_mode()
        del traj

        mtop = md.load_prmtop(topname)
        m_traj = md.load_netcdf(fname, top=mtop)[:1000]

        @Timer()
        def test_rmsd_mdtraj():
            md.rmsd(m_traj, m_traj, 10)

        print("test_rmsd_mdtraj")
        test_rmsd_mdtraj()


if __name__ == "__main__":
    unittest.main()
