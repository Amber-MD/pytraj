import pytest
import pytraj as pt
from pytraj import TrajectoryWriter
from pytraj.testing import aa_eq
from pytraj.utils import tempfolder

# local
from utils import fn

tc5b_trajin = fn('Tc5b.x')
tc5b_top = fn('Tc5b.top')


def test_write_CRYST1():
    traj = pt.datafiles.load_tz2_ortho()[:1]
    print(traj.unitcells)

    with tempfolder():
        fn = "test.pdb"
        traj.save(fn)
        traj2 = pt.load(fn)
        aa_eq(traj.unitcells, traj2.unitcells, decimal=3)


def test_trajectory_writer():
    traj = pt.iterload(tc5b_trajin, tc5b_top)
    with tempfolder():
        pt.write_traj("test_1.pdb", traj[:1], overwrite=True)
        pt.write_traj("test_1.dcd", traj[:1], overwrite=True)

        with TrajectoryWriter("test_1") as trajout:
            trajout.write(traj[0])


def test_write_pdb():
    # TODO: get absolute path so we can use `tempfolder`
    # if not: wrong dir if using TrajectoryIterator
    traj = pt.iterload(tc5b_trajin, tc5b_top)[:]
    TrajectoryWriter()

    # multiple pdb in multiple files, using `save` method in traj
    with tempfolder():
        basename = "test_pdb_files.pdb"
        traj.save(basename, overwrite=True, options="multi")
        for i in range(10):
            fname = basename + "." + str(i + 1)  # cpptraj use `1`
            frame = pt.iterload(fname, traj.top)[0]
            aa_eq(frame.xyz, traj[i].xyz)

    # multiple pdb in multiple files, using `pt.write_traj`
    with tempfolder():
        basename = "test_pdb_files_pt_write_traj.pdb"
        pt.write_traj(basename, traj, overwrite=True, options="multi")
        for i in range(10):
            fname = basename + "." + str(i + 1)  # cpptraj use `1`
            frame = pt.iterload(fname, traj.top)[0]
            aa_eq(frame.xyz, traj[i].xyz)

        with pytest.raises(IOError):
            # write files again, raise if file exists
            pt.write_traj(basename, traj, overwrite=False, options="multi")

    # keepext
    with tempfolder():
        basename = 'test_pdb_files_pt_write_traj.pdb'
        basename2 = 'test_pdb_files_pt_write_traj'
        pt.write_traj(basename, traj, overwrite=True, options="multi keepext")
        for i in range(10):
            fname = "{}.{}.pdb".format(basename2, i+1)
            frame = pt.iterload(fname, traj.top)[0]
            aa_eq(frame.xyz, traj[i].xyz)

        with pytest.raises(IOError):
             pt.write_traj(basename, traj, overwrite=False, options="multi keepext")

    # multiple pdb in SINGLE file
    with tempfolder():
        basename = "test_pdb_files.pdb"
        traj.save(basename, overwrite=True)
        traj2 = pt.load(basename, traj.top)
        aa_eq(traj.xyz, traj2.xyz)

    # multiple pdb in SINGLE file with `model` keyword
    # write to output so we can manually check
    basename = "test_pdb_files_model.pdb"
    with tempfolder():
        traj.save(basename, overwrite=True, options='model')
        traj3 = pt.load(basename, traj.top)
        aa_eq(traj.xyz, traj3.xyz)
