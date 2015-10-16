#!/usr/bin/env python
import pytraj as pt
from memory_profiler import profile

traj = pt.iterload('GAAC3.5000frames.nc', 'GAAC.topo', frame_slice=(0, 5000))

@profile
def calc_pretty(traj, mask='!:WAT', mat_type='cpptraj'):
    return pt.pairwise_rmsd(traj, mask=mask, mat_type=mat_type)

@profile
def calc_not_pretty(traj, mask='!:WAT', mat_type='cpptraj'):
    return pt.pairwise_rmsd(traj(mask=mask), mat_type=mat_type)

calc_pretty(traj, mat_type='full')
calc_not_pretty(traj)
