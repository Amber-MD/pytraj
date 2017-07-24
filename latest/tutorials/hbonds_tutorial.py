
# coding: utf-8

# In[7]:


# load data
# pdb file from: http://www.rcsb.org/pdb/explore/explore.do?pdbId=1l2y
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

import pytraj as pt

traj = pt.iterload('1L2Y.pdb', top='1L2Y.pdb')
traj


# In[10]:


# search hbonds for all residues
# dump data to `dataset` (similiar to python'list and dict)
dataset = pt.hbond(traj)
print(dataset)
# 0: no hbond
# 1: there is hbond


# In[12]:


# get total solute hbonds
print(dataset.total_solute_hbonds())


# In[14]:


# take data for specific hbond
print(dataset.data['ARG16_O-TRP6_NE1-HE1'])

