import pytraj as pt
fn = "../tests/data/Test_NAstruct/adh026.3.pdb"
traj = pt.load(fn, fn)

d = pt.nastruct(traj)
print(d)

# get major groove
print(d.major)

# get minor groove
print(d.minor)
print(d.minor[1].mean(axis=0))

# get inclination
print(d.incl)

# get all supported keys
print(d.keys())
