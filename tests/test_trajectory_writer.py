import unittest
import pytraj as pt
import numpy as np
from pytraj.base import *
from pytraj import io as mdio
from pytraj.testing import aa_eq

farray = pt.load("data/Tc5b.x",
                 "./data/Tc5b.top",
                 frame_indices=list(range(10)))


class TestTrajectoryWriter(unittest.TestCase):

    def test_0(self):
        farray = pt.load("data/Tc5b.x",
                         "./data/Tc5b.top",
                         frame_indices=list(range(10)))
        frame0 = farray[0]
        trajout = TrajectoryWriter()
        trajout.open(filename="./output/test.x",
                     top=farray.top,
                     overwrite=True)
        trajout.write(frame0)

        # add more frames
        for i in range(5, 8):
            trajout.write(farray[i])

        trajout.close()

    def test_1_with_statement(self):
        frame0 = farray[0]
        with TrajectoryWriter(filename="./output/test_trajout_withstatement.x",
                     top=farray.top,
                     overwrite=True) as trajout:
            trajout.write(frame0)

        # reload
        farray2 = Trajectory("./output/test_trajout_withstatement.x",
                             "./data/Tc5b.top")
        frame0_new = farray2[0]

    def test_2(self):
        """test open file writen from test_0"""
        farray = Trajectory()
        farray.top = pt.load_topology('./data/Tc5b.top')
        farray.load("./output/test.x")

    def test_3_write_PDBFILE(self):
        frame0 = farray[0]
        with TrajectoryWriter(filename="./output/test_0.pdb",
                     top=farray.top,
                     overwrite=True) as trajout:
            trajout.write(frame0)

    def test_4(self):
        """test write Trajectory"""
        farray = pt.load("data/Tc5b.x",
                         "./data/Tc5b.top",
                         frame_indices=list(range(10)))
        pt.write_traj("./output/test_write_output.x",
                   farray,
                   top=farray.top,
                   overwrite=True)
        pt.write_traj("./output/test_pdb_1.dummyext",
                   farray[0],
                   top=farray.top,
                   overwrite=True)

        # test 'save'
        farray.save("./output/test_write_output_save_method.x", overwrite=True)

        # reproduce result?
        f0 = mdio.iterload("./output/test_write_output.x", "./data/Tc5b.top")
        f1 = mdio.iterload("./output/test_write_output_save_method.x",
                           "./data/Tc5b.top")
        aa_eq(f0[:, :, :].xyz, f1[:, :, :].xyz)

    def test_5(self):
        farray = Trajectory("./output/test_0.pdb", "./data/Tc5b.top")[0]


if __name__ == "__main__":
    unittest.main()
