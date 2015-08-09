import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        try:
            import mdtraj as md
            has_mdtraj = True
        except ImportError:
            has_mdtraj = False

        if has_mdtraj:
            import mdtraj as md
            traj = traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
            mtop = md.load_prmtop("./data/Tc5b.top")
            print(mtop)
            print (traj.top)

            m_traj = md.Trajectory(traj.xyz, mtop)
            print (m_traj)
            print (traj)

        else:
            print ("does not have mdtraj. skip test")

if __name__ == "__main__":
    unittest.main()
