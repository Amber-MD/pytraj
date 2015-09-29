#!/usr/bin/env python

import pytraj as pt

root_dir = "../../tests/data/nogit/tip3p/"
traj_name = root_dir + "md.nc"
parm_name = root_dir + "tc5bwat.top"

# load to TrajectoryIterator
traj = pt.iterload(traj_name, parm_name, frame_slice=(0, 4000))

state = pt.load_batch(traj, '''
        autoimage
        distance :3 :18
        molsurf @CA
        ''')
state.run()
print(state.data.values)
