#!/usr/bin/env python

import pytraj as pt
from pytraj.testing import Timer

traj = pt.iterload('GAAC3.5000frames.nc', 'GAAC.topo', frame_slice=(0, 1000))
print(traj._estimated_GB)

@Timer()
def iterframe_(traj):
    print('iterframe_')
    for traj in traj:
        pass

@Timer()
def iterchunk_(traj, chunksize):
    print('iterchunk_ : %s' % chunksize)
    for c in pt.iterchunk(traj, chunksize=chunksize):
        pass

iterframe_(traj)
for chunksize in [50, 100, 200, 300, 400, 500, 1000]:
    iterchunk_(traj, chunksize)
