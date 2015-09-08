import pytraj as pt

# use `iterload` to save memory
parm = pt.load_topology("../tests/data/tz2.ortho.parm7")
print(parm.mass)
