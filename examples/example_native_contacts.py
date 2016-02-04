import pytraj as pt

# use iterload for memory saving
traj = pt.iterload("../tests/data/Tc5b.x", "../tests/data/Tc5b.top")

# load reference
frame = pt.load("../tests/data/Tc5b.crd", traj.top)

print(pt.native_contacts(traj, ref=frame))

# get more cpptraj's info (need to use cpptraj keyword)
pt.info("nativecontacts")
