
# coding: utf-8

# In[3]:

import pytraj as pt
traj = pt.iterload('tz2.ortho.nc', 'tz2.ortho.parm7')
#traj = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
traj


# In[7]:

# calculate rmsd to 1st frame, mask = CA
data = pt.rmsd(traj, ref=0, mask='@CA')
print(data)


# In[9]:

# get info about unitcells
traj.unitcells


# In[13]:

# load all coordinates to memory
xyz = traj.xyz[:]
print(xyz[0])


# In[16]:

# calculate bfactor
data = pt.bfactors(traj, '@CA', byres=True)
print(data)


# In[20]:

# plot bfactor
get_ipython().magic('matplotlib inline')
from matplotlib import pyplot as plt
data = data.T
print(data)
plt.plot(data[0], data[1], '--bo')

