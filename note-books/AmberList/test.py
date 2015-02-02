
# coding: utf-8

# In[124]:

# maillist:
# 06/27/2014
from pycpptraj import allactions
from pycpptraj import io as mdio
from pycpptraj.DataSetList import DataSetList
from pycpptraj.DataFileList import DataFileList
from pycpptraj import cast_dataset


# In[125]:

traj = mdio.load("../../tests/data/Ala3/Ala3.crd", "../../tests/data/Ala3/Ala3.top")
dflist = DataFileList()
dslist = DataSetList()


# In[126]:

act = allactions.Action_Matrix()
act.read_input(command="matrix dist @H* out mat.txt", current_top=traj.top, dflist=dflist, dslist=dslist)
act.process(traj.top)


# In[127]:

# distance matrix for first frame
act.do_action(idx=0, current_frame=traj[0])


# In[128]:

dflist.write_all_datafiles()


for idx, atom in enumerate(traj.top['@H*']):
    print idx, atom.name, atom.resnum, traj.top.residuelist[atom.resnum]


# In[137]:

# alternative way
import numpy as np
from pycpptraj.DistRoutines import distance

frame0 = traj[0]
frame0.strip_atoms('!@H*', top=traj.top)

# this is the demo. For large system, use cpptraj or Cython
# for nested loop (which is quite expensive in Python)
# or use enumerate, iterator ...

d0 = cast_dataset(dslist[0], dtype='matrix')
arr1 = np.empty(d0.size)
for i in range(d0.size):
    arr1[i] = d0[i]

print arr1.shape
#arr1.reshape((np.sqrt(d0.size), np.sqrt(d0.size)))
