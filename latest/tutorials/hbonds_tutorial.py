
# coding: utf-8

# In[1]:

# load data
# pdb file from: http://www.rcsb.org/pdb/explore/explore.do?pdbId=1l2y

import pytraj as pt

traj = pt.iterload('1L2Y.pdb', top='1L2Y.pdb')
traj


# In[2]:

# search hbonds for all residues
# dump data to `dataset` (similiar to python'list and dict)
dataset = pt.search_hbonds(traj, dtype='dataset')
print(dataset)
# 0: no hbond
# 1: there is hbond


# In[3]:

# get total solute hbonds
print(dataset['total_solute_hbonds'])


# In[8]:

# take data for specific hbond
print(dataset['ARG16_O-TRP6_NE1-HE1'])


# In[11]:

# take hbond for only ASP9
glu5 = dataset.filter(lambda x: 'ASP9' in x.key)
print(glu5.keys())


# In[18]:

glu5_val = glu5.values
print(glu5_val)

print('total hbonds for GLU5 per frame')
import numpy as np
glu5_total = np.sum(glu5_val, axis=0)
print(glu5_total)


# In[21]:

# plot
get_ipython().magic('matplotlib inline')

from matplotlib import pyplot as plt
plt.plot(glu5_total, 'ro')


# In[22]:

# get usage
help(pt.search_hbonds)

