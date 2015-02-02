from pytraj.base import *
from time import time

traj = Trajectory()
print(traj)
top = Topology("data/Tc5b.top")
traj.top = top
traj.load("data/md1_prod.Tc5b.x")
print(traj)
print(traj.size)
farray = FrameArray()

t0 = time()
for frame in traj:
    farray.append(frame)

print(time() - t0)
print(farray.size)
