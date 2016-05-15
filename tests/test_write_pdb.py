import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import io as mdio
from pytraj import TrajectoryWriter
from glob import glob
from pytraj.testing import aa_eq
from pytraj.utils import tempfolder


class TestWritePDB(unittest.TestCase):

    def test_write_CRYST1(self):
        traj = pt.datafiles.load_tz2_ortho()[:1]
        print(traj.unitcells)

        with tempfolder():
            fn = "test.pdb"
            traj.save(fn)
            traj2 = pt.load(fn)
            aa_eq(traj.unitcells, traj2.unitcells, decimal=3)

    def test_0(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        mdio.write_traj("test_1.pdb", traj[0], top=traj.top, overwrite=True)
        mdio.write_traj("test_1.dcd", traj[0], top=traj.top, overwrite=True)

        with TrajectoryWriter("./output/test_1", overwrite=True) as trajout:
            trajout.write(traj[0])

    def test_1(self):
        # TODO: get absolute path so we can use `tempfolder`
        # if not: wrong dir if using TrajectoryIterator
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")[:]
        trajout = TrajectoryWriter()

        # multiple pdb in multiple files, using `save` method in traj
        with tempfolder():
            basename = "test_pdb_files.pdb"
            traj.save(basename, overwrite=True, options="multi")
            for i in range(10):
                fname = basename + "." + str(i + 1)  # cpptraj use `1`
                frame = mdio.iterload(fname, traj.top)[0]
                aa_eq(frame.xyz, traj[i].xyz)

        # multiple pdb in multiple files, using `mdio.write_traj`
        with tempfolder():
            basename = "test_pdb_files_mdio_write_traj.pdb"
            mdio.write_traj(basename, traj, overwrite=True, options="multi")
            for i in range(10):
                fname = basename + "." + str(i + 1)  # cpptraj use `1`
                frame = pt.iterload(fname, traj.top)[0]
                aa_eq(frame.xyz, traj[i].xyz)

        # multiple pdb in SINGLE file
        with tempfolder():
            basename = "test_pdb_files.pdb"
            traj.save(basename, overwrite=True)
            traj2 = mdio.load(basename, traj.top)
            aa_eq(traj.xyz, traj2.xyz)

        # multiple pdb in SINGLE file with `model` keyword
        # write to output so we can manually check
        basename = "./output/test_pdb_files_model.pdb"
        traj.save(basename, overwrite=True, options='model')
        traj3 = mdio.load(basename, traj.top)
        aa_eq(traj.xyz, traj3.xyz)


if __name__ == "__main__":
    unittest.main()
