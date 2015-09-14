import unittest
from pytraj.base import *
import pytraj.io as mdio
from pytraj.decorators import no_test


class TestPyCpptrajIO(unittest.TestCase):
    def test_save_traj_from_file(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")[:5]
        mdio.write_traj(filename="./output/test_0.binpos",
                        traj=traj,
                        top="./data/Tc5b.top",
                        overwrite=True)

        savedtraj = mdio.iterload("./output/test_0.binpos", traj.top)
        assert savedtraj.n_frames == traj.n_frames

    def test_blindload(self):
        top = mdio.load_topology("./data/Tc5b.top")
        assert isinstance(top, Topology) == True

        traj = mdio.iterload(
            filename="./data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top")

        is_traj = (
            isinstance(traj, TrajectoryIterator) or
            isinstance(traj, Trajectory))
        assert is_traj == True

    def test_ParmFile(self):
        top = mdio.read_parm("./data/Tc5b.top")
        mdio.write_parm("./output/test_io.top", top)
        newtop = mdio.read_parm("./output/test_io.top")
        assert top.n_atoms == newtop.n_atoms

    def test_load_and_save_0(self):
        # need to load to Trajectory to save
        traj = mdio.iterload(
            filename="./data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top")[:]

        indices = list(range(2, 3, 5)) + [3, 7, 9, 8]
        mdio.write_traj(filename="./output/test_io_saved_.x",
                        traj=traj[:],
                        top="./data/Tc5b.top",
                        indices=indices,
                        overwrite=True)

        # check frames
        traj2 = mdio.iterload(
            filename="./output/test_io_saved_.x",
            top="./data/Tc5b.top")

        # about 50% failures
        assert traj2.n_frames == len(indices)

    def test_load_and_save_1(self):
        traj = mdio.iterload(
            filename="./data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top")

        indices = list(range(2, 4)) + [3, 7, 9, 8]
        mdio.write_traj(filename="./output/test_io_saved.pdb",
                        traj=traj,
                        top="./data/Tc5b.top",
                        indices=indices,
                        overwrite=True)

        # check frames
        traj = mdio.iterload(
            filename="./output/test_io_saved.pdb",
            top="./data/Tc5b.top")
        assert traj.n_frames == len(indices)
        assert traj.top.n_atoms == 304


if __name__ == "__main__":
    unittest.main()
