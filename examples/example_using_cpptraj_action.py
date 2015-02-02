import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj import allactions
from pytraj.cast_dataset import cast_dataset

# load traj
farray = mdio.load(top=Topology("../tests/data/Tc5b.top"), 
                       filename='../tests/md1_prod.Tc5b.x',)

# create dataset to hold data
dslist = DataSetList()

# creat Action_Radgyr
act = allactions.Action_Radgyr()

# read input and process input
act.read_input("radgyr @CA", farray.top, dslist=dslist)
act.process(farray.top)

# do calculation for each frame. rad of gyr is appended to dslist
for i, frame in enumerate(farray):
    act.do_action(i, frame)

# get dataset.
d1 = cast_dataset(dslist[0], dtype="general")
print(d1)

# print
print(d1.data[:10])
print(dir(d1))
