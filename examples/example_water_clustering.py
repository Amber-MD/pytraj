import pytraj as pt

# use `iterload` to save memory
# you can use `load`, which is similiar to `mdtraj`
traj = pt.iterload("../tests/data/tz2.ortho.nc",
                   "../tests/data/tz2.ortho.parm7")
# get some info
print(traj)

# get new trajectory for specific waters (Oxygen atom only)
# this command will load all water coordinates to memory
wat_traj = traj[':100-500@O']

# you can use below for lazy-loading
# wat_traj = traj(mask=':100-500@O')
# traj(...) will create iterator, like `range(...)` in python

# iterate every frame and do clustering
for frame in wat_traj:
    xyz = frame.xyz
    # clustering for x-coordniates
    result = pt.clustering_dataset(xyz[:, 0], 'clusters 10')

    # cluster index for each atom
    print(result)
