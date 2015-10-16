#!/usr/bin/env python
import pytraj as pt
import numpy as np
from pytraj.testing import Timer

traj = pt.iterload('GAAC3.5000frames.nc', 'GAAC.topo', frame_slice=(0, 1000))

@Timer()
def iter_speed(xyz, msg):
    print(msg)
    for _ in xyz:
        for _ in xyz:
            for _ in xyz:
                pass

# same
xyz = traj[:].xyz
iter_speed(xyz, 'c_order')

xyz = np.asfortranarray(xyz)
iter_speed(xyz, 'f_order')
