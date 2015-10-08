
# coding: utf-8

# In[27]:

# TODO: validate the data. seems weirds
# make plot inline
get_ipython().magic('matplotlib inline')
from matplotlib import pyplot as plt

# load pytraj
import pytraj as pt

# load traj
traj = pt.load('tz2.ortho.nc', 'tz2.ortho.parm7')
pt.autoimage(traj)
print(traj)


# In[28]:

# do calculation
rdf_data = pt.rdf(traj, bin_spacing=0.2, maximum=20., solvent_mask=':WAT@O', solute_mask=':3')
rdf_data


# In[29]:

# plot
fig = plt.plot(rdf_data[0], rdf_data[1])

