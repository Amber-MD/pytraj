import unittest
import pytraj as pt
import numpy as np
from pytraj.base import *
from pytraj.decorators import no_test
from pytraj.io import write_traj
from pytraj import io as mdio
from pytraj.testing import aa_eq

farray = Trajectory(
    "data/md1_prod.Tc5b.x", "./data/Tc5b.top",
    indices=list(range(10)))


class TestTrajout(unittest.TestCase):
    #@no_test

    def test_0(self):
        farray = Trajectory(
            "data/md1_prod.Tc5b.x", "./data/Tc5b.top",
            indices=list(range(10)))
        frame0 = farray[0]
        trajout = Trajout()
        trajout.open(filename="./output/test.x",
                     top=farray.top,
                     overwrite=True)
        trajout.write(0, frame0, farray.top)

        # add more frames
        for i in range(5, 8):
            trajout.write(i, farray[i], farray.top)

        #assert trajout.is_open() == True
        trajout.close()

    def test_1_with_statement(self):
        frame0 = farray[0]
        with Trajout(filename="./output/test_trajout_withstatement.x",
                     top=farray.top,
                     overwrite=True) as trajout:
            trajout.write(0, frame0, farray.top)

        # reload
        farray2 = Trajectory(
            "./output/test_trajout_withstatement.x", "./data/Tc5b.top")
        frame0_new = farray2[0]

    def test_2(self):
        """test open file writen from test_0"""
        farray = Trajectory()
        farray.top = Topology('./data/Tc5b.top')
        farray.load("./output/test.x")
        #print(farray.size)

        #@no_test
    def test_3_write_PDBFILE(self):
        frame0 = farray[0]
        with Trajout(filename="./output/test_0.pdb",
                     top=farray.top,
                     overwrite=True) as trajout:
            trajout.write(0, frame0, farray.top)

    #@no_test
    def test_4(self):
        """test write Trajectory"""
        farray = Trajectory(
            "data/md1_prod.Tc5b.x", "./data/Tc5b.top",
            indices=list(range(10)))
        write_traj(
            "./output/test_write_output.x", farray, farray.top,
            overwrite=True)
        write_traj("./output/test_pdb_1.dummyext", farray[0], farray.top,
                   overwrite=True)

        # test 'save'
        #print(farray)
        farray.save("./output/test_write_output_save_method.x", overwrite=True)

        # reproduce result?
        f0 = mdio.iterload("./output/test_write_output.x", "./data/Tc5b.top")
        f1 = mdio.iterload(
            "./output/test_write_output_save_method.x", "./data/Tc5b.top")
        aa_eq(f0[:, :, :].xyz, f1[:, :, :].xyz)

    #@no_test
    def test_5(self):
        farray = Trajectory("./output/test_0.pdb", "./data/Tc5b.top")[0]
        #print(farray.n_atoms)


if __name__ == "__main__":
    unittest.main()
