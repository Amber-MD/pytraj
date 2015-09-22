
# coding: utf-8

# PCA analysis with sklearn (change a bit from `mdtraj` tutorial)

# In[1]:

from __future__ import print_function

import pytraj as pt
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from matplotlib.pyplot import rc
font = {'family': 'serif', 'size': '20'}
rc('font', **font)


# In[13]:

# we use `load` method to load all data to memory. This is good for small data size.
# use `pytraj.iterload` for out-of-core traj.

traj = pt.load('data/tz2.nc', 'data/tz2.parm7')
traj


# In[3]:

pca1 = PCA(n_components=2)
pca1


# In[4]:

# superpose to 1st frame
traj.rmsfit(0)


# In[5]:

reduced_cartesian = pca1.fit_transform(traj.xyz.reshape(traj.n_frames, traj.n_atoms * 3))
print(reduced_cartesian.shape)


# In[19]:

plt.figure()
plt.scatter(reduced_cartesian[:, 0], reduced_cartesian[:,1], marker='o', s=30,
        c=range(traj.n_frames))
plt.xlabel('PC1')
plt.ylabel('PC2')
cbar = plt.colorbar()
cbar.set_label('frame #')
plt.savefig('PCA_heart.png')
