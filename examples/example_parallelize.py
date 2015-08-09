import pytraj as pt

# use `iterload` to save memory
traj = pt.iterload("../tests/data/tz2.ortho.nc",
                   "../tests/data/tz2.ortho.parm7")

# use 4 available cores
# calculate radgyr for all atoms
result = pt.pmap(n_cores=4, func=pt.radgyr, traj=traj)
print(result)

# serial version
result = pt.radgyr(traj)
print(result)
