from pytraj import *

traj = FrameArray(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
top = traj.top
frame0 = traj[0]
