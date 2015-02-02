import os
from pytraj.base import *
from pytraj.TrajinList import TrajinList

topname = "./data/Tc5b.top"
trajoutname = "./data/test.x"
refilename = "./data/Tc5b.nat.crd"
trajinname = "./data/md1_prod.Tc5b.x"
toplist = TopologyList()
toplist.add_parm(topname)
toplist.info()

top = toplist[0]

#creat TrajinList instance
trajininput= """
reference Tc5b.nat.crd
"""

argIn = ArgList(trajininput)
trajlist = TrajinList()
trajlist.add_traj("./data/md1_prod.Tc5b.x", argIn, toplist)
trajlist.add_traj("./data/md1_prod.Tc5b.x", argIn, toplist)
trajlist.add_traj("./data/md1_prod.Tc5b.x", argIn, toplist)
trajlist.add_traj("./data/md1_prod.Tc5b.x", argIn, toplist)
print(trajlist.max_frames)
