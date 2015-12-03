import pytraj as pt

traj = pt.load('../tests/data/tz2.nc', '../tests/data/tz2.parm7')

data = pt.pca(traj, mask='@CA', n_vecs=3)
print(pt.pca.__doc__)

print('##################')
print('output')
print(data)

