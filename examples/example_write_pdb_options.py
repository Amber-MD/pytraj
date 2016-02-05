import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj import TrajectoryWriter
from glob import glob
from pytraj.testing import aa_eq
from pytraj.utils import tempfolder


class Test(unittest.TestCase):

    def test_1(self):
        # TODO: get absolute path so we can use `tempfolder`
        # if not: wrong dir if using TrajectoryIterator
        traj = mdio.iterload("../tests/data/Tc5b.x",
                             "../tests/data/Tc5b.top")[:]
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
                frame = mdio.iterload(fname, traj.top)[0]
                aa_eq(frame.xyz, traj[i].xyz)

        # multiple pdb in SINGLE file
        with tempfolder():
            basename = "test_pdb_files.pdb"
            traj.save(basename, overwrite=True)
            traj2 = mdio.load(basename, traj.top)
            aa_eq(traj.xyz, traj2.xyz)

        # multiple pdb in SINGLE file with `optionsl` keyword
        # write to output so we can manually check
        basename = "./output/test_pdb_files_model.pdb"
        traj.save(basename, overwrite=True, options='model')
        traj3 = mdio.load(basename, traj.top)
        aa_eq(traj.xyz, traj3.xyz)


if __name__ == "__main__":
    unittest.main()
