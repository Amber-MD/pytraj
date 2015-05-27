#import unittest
from pytraj import *
from pytraj.utils import Timer
import mdtraj as md


@Timer()
def frame_iter_pytraj():
    traj = io.iterload("md_first_100.nc", "../../tc5bwat.top")
    for frame in traj():
        pass

@Timer()
def chunk_iter_pytraj():
    traj = io.iterload("md_first_100.nc", "../../tc5bwat.top")
    for traj in traj.chunk_iter(chunk=10):
        pass

@Timer()
def chunk_iter_mdtraj():
    for traj in md.iterload("md_first_100.nc", top="tc5bwat.prmtop", chunk=10):
            pass

frame_iter_pytraj()
chunk_iter_pytraj()
chunk_iter_mdtraj()

traj = io.iterload("md_first_100.nc", "../../tc5bwat.top")
print (traj)

# result:
#    0.095
#    0.180
#    1.852
#    TrajectoryIterator instance with 100 frames, 31712 atoms/frame
