import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj import Trajout

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj import set_world_silent
        set_world_silent(False)
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        mdio.writetraj("test_1.pdb", traj[0], top=traj.top, fmt='CHARMMDCD', overwrite=True)
        mdio.writetraj("test_1.dcd", traj[0], top=traj.top, fmt='CHARMMDCD', overwrite=True)

        with Trajout("./output/test_1", fmt="PDBFILE", overwrite=True) as trajout:
            trajout.writeframe(frame=traj[0], top=traj.top)

    def test_1(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        trajout = Trajout()
        print(trajout.formats)
         
if __name__ == "__main__":
    unittest.main()
