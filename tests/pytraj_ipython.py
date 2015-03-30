from pytraj import *
from pytraj.common_actions import *
from pytraj.datasets import DataSet_Coords_TRJ

traj = FrameArray(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
traj2 = io.load('./data/DPDP.nc', "./data/DPDP.parm7")
top = traj.top
frame0 = traj[0]

act = adict['matrix']
dslist = DataSetList()
act("byres @CA", traj, dslist=dslist)

coords_traj = DataSet_Coords_TRJ()
coords_traj.top = traj2.top
coords_traj.add_trajin(traj2)
