import unittest  # for travis (you don't need to add this)

# since MD simulations produce huge amount of data (hundreds of GB, few TB
# or enven tens of TB)
# so we always try to use TrajectoryIterator to save memory
# coordinates of specific frames are only loaded if needed

import pytraj as pt
traj = pt.iterload('../tests/data/tz2.ortho.nc',
                   '../tests/data/tz2.ortho.parm7')

# calculate molsurf for frame 0, 2, 4, 6
# start=0, stop=8, stride=2 (just like python's range(0, 8, 2))
print(pt.molsurf(traj(0, 8, 2), '@CA'))

# if we want to load all frames at once, use [] indexing style
# (just like indexing a list/array in python)
print(pt.molsurf(traj[0:8:2], '@CA'))

# we even perform `autoimage`
print(pt.radgyr(traj(0, 8, 2, autoimage=True), '@CA', nomax=True))
