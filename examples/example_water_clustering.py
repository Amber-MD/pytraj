import pytraj as pt

# use `iterload` to save memory
# you can use `load`, which is similiar to `mdtraj`
traj = pt.iterload("../tests/data/tz2.ortho.nc",
                   "../tests/data/tz2.ortho.parm7")
# get some info
print(traj)

# get new trajectory for specific waters (Oxygen atom only)
wat_traj = traj[':100-500@O']

# iterate every frame and do clustering
for frame in wat_traj:
    xyz = frame.xyz
    # clustering for x-coordniates
    result = pt.clustering_dataset(xyz[:, 0], 'clusters 10')

    # cluster index for each atom
    print(result)
