#!/usr/bin/env python
import pytraj as pt
from memory_profiler import profile

traj = pt.iterload('GAAC3.5000frames.nc', 'GAAC.topo', frame_slice=(0, 5000))

@profile
def calc_(traj=traj):
    frame = pt.mean_structure(traj)
    print(frame)

calc_()

# output:
# Line #    Mem usage    Increment   Line Contents
# ================================================
#      7     67.9 MiB      0.0 MiB   @profile
#      8                             def calc_(traj=traj):
#      9    111.8 MiB     43.9 MiB       frame = pt.mean_structure(traj)
#     10    111.8 MiB      0.0 MiB       print(frame)
