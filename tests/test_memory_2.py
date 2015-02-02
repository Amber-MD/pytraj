import numpy as np
#from memory_profiler import profile
#from line_profiler import profile
from pytraj.base import *

TRAJ = TrajReadOnly()
TRAJ.top = Topology("./data/Tc5b.top")
TRAJ.load("./data/md1_prod.Tc5b.x")

#@profile
def calc_pairwise_rmsd():
    farray = FrameArray()
    farray.top = TRAJ.top
    #
    for frame in TRAJ:
        frame.strip_atoms("!@CA", TRAJ.top.copy())
        farray.append(frame)
    
    size = farray.size
    arr = np.empty(shape=(size, size))
    #
    for i, frame_i in enumerate(farray):
        for j, frame_j in enumerate(farray):
            arr[i, j] = frame_i.rmsd(frame_j)

if __name__ == "__main__":
    calc_pairwise_rmsd()
