
# coding: utf-8

# <span style="color:blue; font-family:Georgia; font-size:3em;"> **Calculate pairwise RMSD for trajectory and plot** </span>

# In[1]:

import pytraj as pt
from pytraj.plot import plot_matrix


# In[2]:

traj = pt.iterload("../tests/data/md1_prod.Tc5b.x", "../tests/data/Tc5b.top")


# In[3]:

mat = pt.pairwise_rmsd(traj, '@CA')
mat


# In[4]:

get_ipython().magic('matplotlib inline')
get_ipython().magic("config InlineBackend.figure_format = 'retina'")
import matplotlib
matplotlib.rcParams['savefig.dpi'] = 2 * matplotlib.rcParams['savefig.dpi']

fig, asp, axi, = plot_matrix(mat)
fig.colorbar(axi, ax=asp)

