import unittest
from pytraj.base import *
from pytraj import io as mdio

class Test(unittest.TestCase):
    def test_0(self):
        # TODO : not right yet.
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        # this command actually write Amber restart file
        #mdio.writetraj("test_1.pdb", traj[0], top=traj.top, fmt='CHARMM')
        print(dir(Trajout))
        print(Trajout.writeframe)
        #with Trajout("test_1.pdb", fmt='PDBFILE', overwrite=True) as trajout:
        #with Trajout("test_1.x", fmt='MOL2FILE', overwrite=True) as trajout:
        with Trajout("./output/test_1", fmt="PDBFILE", overwrite=True) as trajout:
            trajout.writeframe(frame=traj[0], top=traj.top)

    def test_1(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        trajout = Trajout()
        print(trajout.formats)
         
if __name__ == "__main__":
    unittest.main()
