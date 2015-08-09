import pytraj as pt

# use `iterload` to save memory
traj = pt.iterload("../tests/data/tz2.ortho.nc",
                   "../tests/data/tz2.ortho.parm7")
print(traj)

# load water-stripped trajectory to memory
t0 = traj['!:WAT']
t0.top.set_nobox()
print(t0)
