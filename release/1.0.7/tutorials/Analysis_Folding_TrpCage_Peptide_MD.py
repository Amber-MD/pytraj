
# coding: utf-8

# ## Trajectory analysis  for MD Simulation from Folding Trp-Cage Peptide
# 
#     Reproduce some analysis in http://ambermd.org/tutorials/basic/tutorial3/section6.htm

# ## Download data
# 
# ```
# url: http://ambermd.org/tutorials/basic/tutorial3/section6.htm
# ```
# 
# ### Require
# - pytraj/cpptraj
# - nglview (https://github.com/arose/nglview)
#     ```bash
#           # Linux users
#           conda install -c ambermd nglview
#           
#           # for Mac user, please see the instruction from above website
#     ```
#     
# - pandas
# 
#     ```bash
#          conda install pandas
#     ```

# In[1]:

# uncomment below to download, untar files.

#!wget http://ambermd.org/tutorials/basic/tutorial3/files/production.tar.gz
#!wget http://ambermd.org/tutorials/basic/tutorial3/files/TC5b.prmtop
#!mv TC5b.prmtop production/
#!tar -xf ./production.tar.gz
#! wget http://ambermd.org/tutorials/basic/tutorial3/files/lowest_energy_struct.pdb.2455 production/
# !mv lowest_energy_struct.pdb.2455 production/
get_ipython().system('ls production/')


# ## Load trajectory files by pytraj

# In[2]:

import pytraj as pt
traj = pt.load('production/equil*.mdcrd.gz', top='production/TC5b.prmtop')
traj


# ## Visualization with nglview

# In[3]:

# note: you need to run the notebook
import nglview as nv

view = nv.show_pytraj(traj)
view


# In[4]:

# reset representation
view.representations = []
view.parameters = {'theme': 'dark'}
view.add_representation('licorice', selection='not hydrogen')
view.add_representation('cartoon')


# In[5]:

# jump to specific frame
view.frame = 2000


# In[6]:

# if you're using notebook, you should see somethinge similiar to below image
from IPython.display import Image
Image('nglview_trpcage_amber_tutorial_3.png', width=600)


# ## Plot RMSD vs time to lowest structure (taken from tutorial)

# In[7]:

# load reference (we don't need to use a topology file since this is pdb)
lowest_energy_pdb = pt.load('production/lowest_energy_struct.pdb.2455')
lowest_energy_pdb


# In[8]:

rmsd_to_lowest_pdb = pt.rmsd(traj, ref=lowest_energy_pdb, mask='@N,CA,C')
rmsd_to_lowest_pdb


# In[9]:

get_ipython().magic('matplotlib inline')

from matplotlib import pyplot as plt

plt.plot(rmsd_to_lowest_pdb)
plt.xlabel('Time (ps)')
plt.ylabel('RMSD (angstrom)')


# ## Compute dihedral angle for TRP6

# In[10]:


# store data in pandas's DataFrame
# 
df_trp6_dihedrals = pt.multidihedral(traj, resrange='6', dtype='dataframe')
df_trp6_dihedrals.head(6)


# In[11]:

plt.plot(df_trp6_dihedrals['phi_6'], 'bo', markersize=1.)


# # Hbond

# In[12]:

hbond_data = pt.hbond(traj)
hbond_data


# In[13]:

## Get donor and acceptor mask, only show first 10
hbond_data.donor_aceptor[:10]


# ## plot total hbond number vs rmsd

# In[19]:

# plot total hbond number vs rmsd
# color by frame number
n_hbonds_per_frame = hbond_data.total_solute_hbonds()
fig = plt.scatter(n_hbonds_per_frame, rmsd_to_lowest_pdb, 
                  marker='o', c=range(traj.n_frames), alpha=0.8, cmap='Spectral')
plt.colorbar()


# In[24]:

# get distance and angle mask for each hbond
h_dist_mask, h_angle_mask = hbond_data.get_amber_mask()
print(h_dist_mask)
print(h_angle_mask)


# In[30]:

# compute hbond distance
# use dtype='dataframe' to dump distance data to pandas' DataFrame
h_dist = pt.distance(traj, h_dist_mask, dtype='dataframe')

# need to update mask
h_dist.columns = h_dist_mask


# In[32]:

# just want to show some data
h_dist[[':1@OD1 :2@H', ':8@O :11@H', ':16@O :16@HE']].head(5)


# In[34]:

# stats
h_dist.describe()


# In[37]:

# plot
h_dist[':1@OD1 :2@H'].hist()


# In[38]:

h_dist[':1@OD1 :2@H'].describe()

