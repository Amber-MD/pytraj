# require ParmEd
# pip install git+git://github.com/ParmEd/ParmEd

import pytraj as pt

# load as trajectory
traj = pt.load(
    "http://ambermd.org/tutorials/basic/tutorial4/files/sustiva.mol2")
print(traj)
