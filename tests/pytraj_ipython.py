from pytraj import *
from pytraj.common_actions import *

traj = FrameArray(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
traj2 = io.load('./data/DPDP.nc', "./data/DPDP.parm7")
top = traj.top
frame0 = traj[0]
