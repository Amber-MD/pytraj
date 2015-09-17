#!/usr/bin/env python

import pytraj.io as io
from pytraj import set_world_silent
from pytraj import dihedral_analysis as da
from pytraj.TrajectoryIterator import TrajectoryIterator

set_world_silent(False)

root_name = "../../md"
traj = TrajectoryIterator()
traj.top = io.load_topology("../../prmtop")
print(traj)

for i in range(1, 22):
    fname = root_name + str(i) + ".nc"
    traj.load(fname)
print(traj)

iter_options = {'stop': -1, 'stride': 1}
dslist = da.calc_multidihedral(traj(**iter_options), "chin alpha gamma",
                               top=traj.top)

io.to_pickle(dslist.to_dict(), "chin_alpha_gamma.pk")
#print (dslist.to_dict())
