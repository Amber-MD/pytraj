
# coding: utf-8

# ### PCA analysis with sklearn (adapted from `mdtraj` and cpptraj tutorial)

# In[1]:

import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

get_ipython().magic('matplotlib inline')
from __future__ import print_function
get_ipython().magic("config InlineBackend.figure_format = 'retina'  # high resolution")
import matplotlib
matplotlib.rcParams['savefig.dpi'] = 2 * matplotlib.rcParams['savefig.dpi'] # larger image

import pytraj as pt
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from pytraj.plot import show_config


# In[2]:

# we use `load` method to load all data to memory. This is good for small data size.
# use `pytraj.iterload` for out-of-core traj.

traj = pt.load('tz2.nc', 'tz2.parm7')
traj


# In[3]:

pca = PCA(n_components=2)
pca


# In[4]:

# superpose to 1st frame
pt.superpose(traj, ref=0, mask='!@H=')


# In[5]:

# create average structure

avg = pt.mean_structure(traj)
avg


# In[6]:

# superpose all structures to average frame
pt.superpose(traj, ref=avg, mask='!@H=')


# In[ ]:

# perform PCA calculation and get transformed coords
# we need to reshape 3D traj.xyz array to 2D to make sklearn happy
# make a new traj by stripping all H atoms
traj_new = traj['!@H=']
xyz_2d = traj_new.xyz.reshape(traj_new.n_frames, traj_new.n_atoms * 3)
print(xyz_2d.shape) # (n_frames, n_dimensions)


# In[ ]:

reduced_cartesian = pca.fit_transform(xyz_2d)
print(reduced_cartesian.shape) # (n_frames, n_dimensions)


# In[ ]:

# ignore warning
import warnings
warnings.filterwarnings('ignore')

plt.figure()
plt.scatter(reduced_cartesian[:, 0], reduced_cartesian[:,1], marker='o', c=range(traj_new.n_frames), alpha=0.5)
plt.xlabel('PC1')
plt.ylabel('PC2')
cbar = plt.colorbar()
cbar.set_label('frame #')


# ### Compare to cpptraj data
# 
# **note**: stop here if you do not care (a bit compilicated code)

# In[ ]:

# cpptraj

# copy from Amber15 manual (page 619)
command = '''
# Step one. Generate average structure.
# RMS-Fit to first frame to remove global translation/rotation.
parm tz2.parm7
trajin tz2.nc
rms first !@H=
average crdset AVG
run
# Step two. RMS-Fit to average structure. Calculate covariance matrix.
# Save the fit coordinates.
rms ref AVG !@H=
matrix covar name MyMatrix !@H=
createcrd CRD1
run
# Step three. Diagonalize matrix.
runanalysis diagmatrix MyMatrix vecs 2 name MyEvecs
# Step four. Project saved fit coordinates along eigenvectors 1 and 2
crdaction CRD1 projection evecs MyEvecs !@H= out project.dat beg 1 end 2
'''


# In[ ]:

state = pt.datafiles.load_cpptraj_state


# In[ ]:

state = pt.datafiles.load_cpptraj_state(command)
# tell 'run' to perform all calculation
state.run()


# In[ ]:

# get data
state.data


# In[ ]:

print([dset.key for dset in state.data])
print(state.data['MyMatrix'].values.shape)


# In[ ]:

# reduced_cartesian corresponds to dataset with names of 'Mode1', 'Mode2'
mode_0, mode_1 = -state.data['Mode1'].values, -state.data['Mode2'].values
# mode_0, mode_1 = state.data['Mode1'].values, state.data['Mode2'].values

# plot: cpptraj
fig = plt.figure()
ax_0 = fig.add_subplot(211)
ax_0.scatter(mode_0, mode_1, marker='o', c=range(traj.n_frames), alpha=0.5)
ax_0.set_xlabel('PC1')
ax_0.set_ylabel('PC2')
ax_0.set_xlim([-60, 80])
ax_0.set_ylim([-40, 40])
ax_0.set_title('cpptraj-flip')
ax_0.set_yticks([-40, -20, 0, 20, 40])

# plot: sklearn
fig = plt.figure()
ax_1 = fig.add_subplot(212)
ax_1.scatter(reduced_cartesian[:, 0], reduced_cartesian[:,1], marker='o', c=range(traj.n_frames), alpha=0.5)
ax_1.set_xlabel('PC1')
ax_1.set_ylabel('PC2')
ax_1.set_yticks([-40, -20, 0, 20, 40])
ax_1.set_title('sklearn')


# In[ ]:

print('sklearn \n')
print(reduced_cartesian[:, 0], reduced_cartesian[:, 1])

print('\ncpptraj\n') # flip sign, why?
print(mode_0, mode_1)

