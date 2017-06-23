
# coding: utf-8

# ### load trajectory

# In[1]:

import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)
# load pytraj
import pytraj as pt

# load trajectory to memory
traj = pt.load('tz2.nc', 'tz2.parm7')


# ### Get the data

# In[2]:

# compute pca
data = pt.pca(traj, mask='!@H=', n_vecs=2)

print('projection values of each frame to first mode = {} \n'.format(data[0][0]))
print('projection values of each frame to second mode = {} \n'.format(data[0][1]))
print('eigvenvalues of first two modes', data[1][0])
print("")
print('eigvenvectors of first two modes: \n', data[1][1])


# ### plot

# In[3]:

get_ipython().magic('matplotlib inline')
get_ipython().magic("config InlineBackend.figure_format = 'retina'  # high resolution")
import matplotlib

projection_data = data[0]
from matplotlib import pyplot as plt

plt.scatter(projection_data[0], projection_data[1], marker='o', c=range(traj.n_frames), alpha=0.5)
plt.xlabel('PC1')
plt.ylabel('PC2')
cbar = plt.colorbar()
cbar.set_label('frame #')

