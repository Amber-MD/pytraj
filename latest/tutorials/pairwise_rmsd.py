
# coding: utf-8

# # Calculate pairwise RMSD for trajectory and plot

# In[1]:

import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

import pytraj as pt
from pytraj.plot import plot_matrix


# In[2]:

traj = pt.iterload("tz2.nc", "tz2.parm7")


# In[3]:

mat = pt.pairwise_rmsd(traj, '@CA')
mat


# In[4]:

get_ipython().magic('matplotlib inline')
get_ipython().magic("config InlineBackend.figure_format = 'retina'")
import matplotlib

fig, asp, axi, = plot_matrix(mat)
fig.colorbar(axi, ax=asp)

# from matplotlib import pyplot as plt
# plt.savefig('plot_pairwise_rmsd.png')

