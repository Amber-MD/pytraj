import pytraj as pt

# use iterload for memory saving
# you can use `load` to load all data to mem too

traj0 = pt.iterload(['../tests/data/Test_RemdTraj/rem.nc.000',
                     '../tests/data/Test_RemdTraj/rem.nc.001'],
                    top='../tests/data/Test_RemdTraj/ala2.99sb.mbondi2.parm7')

print(traj0)

# you can use `*` to load all files too
traj1 = pt.iterload('../tests/data/Test_RemdTraj/rem.nc.*',
                    top='../tests/data/Test_RemdTraj/ala2.99sb.mbondi2.parm7')

print(traj1)

# just want to calculate radgyr for few frames?
# frames: 1, 3, 6, 8
print(pt.radgyr(traj1[[1, 3, 6, 8]]))

# frame 1 to 100, skip every 2 frames
print(pt.radgyr(traj1(1, 100, 2)))
