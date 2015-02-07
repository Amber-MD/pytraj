import unittest
from pytraj.base import *
import pytraj.io as mdio
from pytraj.decorators import no_test

class TestPyCpptrajIO(unittest.TestCase):
    @no_test
    def test_save_traj_from_file(self):
        print("test_save_traj_from_file")
        Trajout().help()
        traj = mdio.load("./data/md1_prod.Tc5b.x",
                        "./data/Tc5b.top")[:5]
        print(traj.size)
        mdio.writetraj(filename="./output/test_0.binpos", 
                       traj=traj, 
                       top="./data/Tc5b.top",
                       fmt="BINPOS",
                       overwrite=True)

        savedtraj = mdio.load("./output/test_0.binpos", traj.top)
        print("test_0.binpos size = ", savedtraj.size)
        print(traj.size)
        assert savedtraj.size == traj.size

    def test_blindload(self):
        print("test_blindload")
        top = mdio.load("./data/Tc5b.top")
        assert isinstance(top, Topology) == True

        traj = mdio.load(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
        print(traj)

        # TODO : add more classes here
        is_traj = (isinstance(traj, TrajReadOnly) or isinstance(traj, FrameArray))
        assert is_traj == True

    def test_ParmFile(self):
        print("test_ParmFile")
        top = mdio.readparm("./data/Tc5b.top")
        mdio.writeparm("./output/test_io.top", top)
        newtop = mdio.readparm("./output/test_io.top")
        assert top.n_atoms == newtop.n_atoms

    def test_load_and_save_0(self):
        print("test_load_and_save_0")
        # need to load to FrameArray to save
        traj = mdio.load(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")[:]
        trajiter = mdio.iterload(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top", chunk=100)
        print(trajiter)
        print(traj.size)

        indices = list(range(2, 3, 5)) + [3, 7, 9, 8]
        print(type(indices))
        print(indices)
        mdio.writetraj(filename="./output/test_io_saved_.x", 
                       traj=traj[:], 
                       top="./data/Tc5b.top",
                       indices=indices,
                       overwrite=True)

        # check frames
        traj2 = mdio.load(filename="./output/test_io_saved_.x", top="./data/Tc5b.top")
        print("test_load_and_save_0")
        print(traj2.size)
        print(traj2.is_empty())
        print(len(indices))

        # about 50% failures
        # TODO : why?
        assert traj2.size == len(indices) 

    def test_load_and_save_1(self):
        print("test_load_and_save_1")
        traj = mdio.load(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
        trajiter = mdio.iterload(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top", chunk=9)
        print(trajiter)

        indices = list(range(2, 4)) + [3, 7, 9, 8]
        mdio.writetraj(filename="./output/test_io_saved.pdb", 
                       traj=traj, 
                       top="./data/Tc5b.top",
                       fmt='pdbfile', 
                       indices=indices,
                       overwrite=True)

        # check frames
        traj = mdio.load(filename="./output/test_io_saved.pdb", top="./data/Tc5b.top")
        assert traj.size == len(indices) 
        assert traj.top.n_atoms == 304

if __name__ == "__main__":
    unittest.main()
