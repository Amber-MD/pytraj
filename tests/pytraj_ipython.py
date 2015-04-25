from pytraj import *
from pytraj.common_actions import *
from pytraj.datasets import DataSet_Coords_TRJ
from pytraj.TrajinList import TrajinList

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

import mdtraj as md
import mdtraj.testing
m_traj = md.load("./data/DPDP.nc", top="./data/DPDP.parm7")

trajlist = TrajinList()
trajlist.add_traj("./data/md1_prod.Tc5b.x", top, "1 3")
trajlist.add_traj("./data/md1_prod.Tc5b.x", top, "4 9 2")
trajlist.add_traj("./data/md1_prod.Tc5b.x", top, "5 7")
trajlist.add_traj("./data/md1_prod.Tc5b.x", top, "1 last")
