import pytraj as pt

from pytraj.testing import aa_eq, tempfolder
from pytraj import *

# local
from utils import fn
farray = pt.load(fn("Tc5b.x"),
                 fn("Tc5b.top"),
                 frame_indices=list(range(10)))


def test_trajectory_writer_open_close():
    farray = pt.load(fn("Tc5b.x"),
                     fn("Tc5b.top"),
                     frame_indices=list(range(10)))
    frame0 = farray[0]
    with tempfolder():
        trajout = TrajectoryWriter()
        trajout.open(filename="test.x",
                     top=farray.top,
                     overwrite=True)
        trajout.write(frame0)

        # add more frames
        for i in range(5, 8):
            trajout.write(farray[i])

        trajout.close()

        farray = Trajectory()
        farray.top = pt.load_topology(fn('Tc5b.top'))
        farray.load("test.x")

def test_trajectory_writer__with_statement():
    frame0 = farray[0]
    with tempfolder():
        with TrajectoryWriter(filename="test_trajout_withstatement.x",
                     top=farray.top,
                     overwrite=True) as trajout:
            trajout.write(frame0)
        # reload
        farray2 = Trajectory("test_trajout_withstatement.x", fn('Tc5b.top'))
        farray2[0]

def test_trajectory_writer_write_PDBFILE():
    frame0 = farray[0]
    with TrajectoryWriter(filename="./output/test_0.pdb",
                 top=farray.top,
                 overwrite=True) as trajout:
        trajout.write(frame0)

def test_trajectory_writer_write_Trajectory():
    """test write Trajectory"""
    farray = pt.load(fn("Tc5b.x"),
                     fn("Tc5b.top"),
                     frame_indices=list(range(10)))
    with tempfolder():
        pt.write_traj("test_write_output.x",
                   farray,
                   top=farray.top,
                   overwrite=True)
        pt.write_traj("test_pdb_1.dummyext",
                   farray[0],
                   top=farray.top,
                   overwrite=True)

        # test 'save'
        farray.save("test_write_output_save_method.x", overwrite=True)

        # reproduce result?
        f0 = pt.iterload("test_write_output.x", fn('Tc5b.top'))
        f1 = pt.iterload("test_write_output_save_method.x", fn('Tc5b.top'))
        aa_eq(f0[:, :, :].xyz, f1[:, :, :].xyz)
