#!/usr/bin/env python

import pytraj as pt
from pytraj.testing import Timer

traj = pt.iterload('GAAC3.5000frames.nc', 'GAAC.topo', frame_slice=(0, 1000))

@Timer()
def pick_topology_(top):
    pt.to_pickle(top, 'test.pk')

for _ in range(5):
    pick_topology_(traj.top)
