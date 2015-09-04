#!/usr/bin/env python
import pytraj as pt
from memory_profiler import profile

traj0 = pt.iterload('GAAC3.5000frames.nc', 'GAAC.topo', frame_slice=(0, 50))
traj1 = pt.iterload('GAAC3.5000frames.nc', 'GAAC.topo', frame_slice=(100, 150))

frame_indices = range(50)


def calc(traj0, traj1, mask='@P'):
    import numpy as np
    arr = np.empty((traj0.n_frames, traj1.n_frames), dtype='f8')

    for idx, frame in enumerate(traj1):
        arr[idx] = pt.rmsd(traj0, ref=frame, mask=mask)
    return arr


print(calc(traj0, traj1))
