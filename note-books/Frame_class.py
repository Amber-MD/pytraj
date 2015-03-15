from __future__ import print_function
# coding: utf-8

## Frame object

# Just like cpptraj, Frame object is work-horse of pytraj. 

# In[53]:

from pytraj import io as mdio
from pytraj import Frame


# In[54]:

# loading traj file to FrameArray object
traj = mdio.load("../tests/data/md1_prod.Tc5b.x", "../tests/data/Tc5b.top")


# In[55]:

# get frame object
frame0 = traj[0]

# frame0 behaves like 2D array with shape of (n_atoms, 3)
print(frame0.shape)


# In[56]:

# fancy indexing
# whole coords
frame0[:][0]


# In[57]:

# or
frame0[0, :]


# In[58]:

# use numpy array as memory view for Frame object
import numpy as np

arr0 = np.asarray(frame0[:])

# update arr0 will update frame coords
print(frame0[0, 0])
arr0[0, 0] = 1000.
print(frame0[0, 0])


# In[59]:

# extracting Frame coords with given mask. 
# 1st way

print(frame0[traj.top("@CA")])


# In[60]:

# 2nd way
# need to set Topology for frame to use AtomMask
frame0.set_top(traj.top)
frame0["@CA"]


# In[61]:

# do basic math with Frame object (you can use numpy memory (as demonstated before))

print(frame0[12, :])
frame0 += frame0
print(frame0[12, :])


# In[62]:

# calculate rmsd between two frames (2-th frame and 9-th frame in traj object)

print(traj[2].rmsd(traj[9]))


# In[63]:

# methods / properties
print(dir(frame0))


### Perform cpptraj Action on Frame object

# In[64]:

# this code shows how we're able control workflow of cpptraj.
from pytraj import allactions
radgyr = allactions.Action_Radgyr()


# In[65]:

# import DataSetList object to store radius of gyration data
# there is much shorter way to do this but this is for demonstration
from pytraj.DataSetList import DataSetList

# store data file for writing output
from pytraj.DataFileList import DataFileList


# In[66]:

dsetlist = DataSetList()
dflist = DataFileList()
# perform action on traj

# calculate radgyr using CA, 
radgyr.read_input("radgyr @CA out test.out", current_top=traj.top, dslist=dsetlist, dflist=dflist)


# In[67]:

# process Topology if needed
radgyr.process(traj.top)


# In[68]:

# start looping all Frame objects


# In[69]:

for idx, frame in enumerate(traj):
    radgyr.do_action(idx, frame)


# In[70]:

# it's time to get the data
import pytraj as pyc

# currently we need to cast dataset since cpptraj has several kinds
d0 = pyc.cast_dataset(dsetlist[0])
print(d0[:])


# In[71]:

# saving to file
dflist.write_all_datafiles()


# In[72]:

# make sure we already save it
get_ipython().system('ls test.out')


# In[73]:

get_ipython().system('head test.out')


# In[73]:



