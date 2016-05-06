import pytraj as pt
from pytraj.testing import get_fn, aa_eq

fn, tn = get_fn('tz2_dry')

traj0 = pt.load(fn, tn)
traj1 = pt.load(fn, tn)
traj2 = pt.load(fn, tn)
print(traj0.xyz[0, 0], traj1.xyz[0, 0], traj2.xyz[0, 0])

# Use ActionList for traj0
actlist = pt.ActionList(top=traj0.top)
actlist.add("translate", "x 1.2")
actlist.add("center", "origin")
actlist.add("rotate", "x 45.")

for frame in traj0:
    actlist.compute(frame)

# use transformation
itertraj = pt.transform(traj1, by=['translate x 1.2', 'center origin', 'rotate x 45.'])
for frame in itertraj:
    pass

aa_eq(traj0.xyz, traj1.xyz)
print(traj0.xyz[0, 0], traj1.xyz[0, 0], traj2.xyz[0, 0])

# use API
traj2.translate('x 1.2')
traj2.center('origin')
traj2.rotate('x 45.')

aa_eq(traj0.xyz, traj2.xyz)
print(traj0.xyz[0, 0], traj1.xyz[0, 0], traj2.xyz[0, 0])
