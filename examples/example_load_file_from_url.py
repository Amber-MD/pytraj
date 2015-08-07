# require ParmEd
# pip install git+git://github.com/ParmEd/ParmEd

import pytraj as pt

# load as trajectory
traj = pt.load(
    "http://ambermd.org/tutorials/advanced/tutorial1/files/polyAT_edit.pdb")
print(traj)

traj = pt.load(
    "http://ambermd.org/tutorials/basic/tutorial4/files/sustiva.mol2")
print(traj)

# load as Toplogy file
top = pt.load_topology(
    "http://ambermd.org/tutorials/advanced/tutorial1/files/polyAT_edit.pdb")
print(top)

# mol2 file
top = pt.load_topology(
    "http://ambermd.org/tutorials/basic/tutorial4/files/sustiva.mol2")
print(top)
