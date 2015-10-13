#!/usr/bin/env python

import pytraj as pt

root_dir = "../../tests/data/"
traj_name = root_dir + "tz2.nc"
parm_name = root_dir + "tz2.parm7"

# load to TrajectoryIterator
traj = pt.iterload(traj_name, parm_name, frame_slice=(0, 4000))

state = pt.load_batch(traj, '''
        autoimage
        distance :3 :10
        molsurf @CA
        ''')

state.run()
print([d for d in state.data])
