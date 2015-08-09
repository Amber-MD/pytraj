import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj import Trajout
from glob import glob
from pytraj.testing import aa_eq
from pytraj.utils import goto_temp_folder


class Test(unittest.TestCase):
    def test_0(self):
        from pytraj import set_world_silent
        set_world_silent(False)
        traj = mdio.iterload("../tests/data/md1_prod.Tc5b.x",
                             "../tests/data/Tc5b.top")
        mdio.write_traj("test_1.pdb", traj[0],
                        top=traj.top,
                        format='CHARMMDCD',
                        overwrite=True)
        mdio.write_traj("test_1.dcd", traj[0],
                        top=traj.top,
                        format='CHARMMDCD',
                        overwrite=True)

        with Trajout("./output/test_1",
                     format="PDBFILE",
                     overwrite=True) as trajout:
            trajout.write(frame=traj[0], top=traj.top)

    def test_1(self):
        # TODO: get absolute path so we can use `goto_temp_folder`
        # if not: wrong dir if using TrajectoryIterator
        traj = mdio.iterload("../tests/data/md1_prod.Tc5b.x",
                             "../tests/data/Tc5b.top")[:]
        trajout = Trajout()
        print(trajout.formats)

        # multiple pdb in multiple files, using `save` method in traj
        with goto_temp_folder():
            basename = "test_pdb_files.pdb"
            traj.save(basename, overwrite=True, mode="multi")
            for i in range(10):
                fname = basename + "." + str(i + 1)  # cpptraj use `1`
                print(fname)
                frame = mdio.load(fname, traj.top)[0]
                aa_eq(frame.xyz, traj[i].xyz)

        # multiple pdb in multiple files, using `mdio.write_traj`
        with goto_temp_folder():
            basename = "test_pdb_files_mdio_write_traj.pdb"
            mdio.write_traj(basename, traj, overwrite=True, mode="multi")
            for i in range(10):
                fname = basename + "." + str(i + 1)  # cpptraj use `1`
                print(fname)
                frame = mdio.load(fname, traj.top)[0]
                aa_eq(frame.xyz, traj[i].xyz)

        # multiple pdb in SINGLE file
        with goto_temp_folder():
            basename = "test_pdb_files.pdb"
            traj.save(basename, overwrite=True)
            traj2 = mdio.load(basename, traj.top)
            aa_eq(traj.xyz, traj2.xyz)

        # multiple pdb in SINGLE file with `model` keyword
        # write to output so we can manually check
        basename = "./output/test_pdb_files_model.pdb"
        traj.save(basename, overwrite=True, mode='model')
        traj3 = mdio.load(basename, traj.top)
        aa_eq(traj.xyz, traj3.xyz)


if __name__ == "__main__":
    unittest.main()
