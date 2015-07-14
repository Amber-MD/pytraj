import pytraj as pt

# use `iterload` to save memory
parm = pt.Topology("../tests/data/tz2.ortho.parm7")
print (parm.mass)
