
# coding: utf-8

# # require
# 
# - pytraj/cpptraj: http://amber-md.github.io/pytraj/latest/installation.html#from-source-code
# 
#         git clone https://github.com/Amber-MD/pytraj
#         cd pytraj
#         python setup.py install
#         
# - libsander (http://ambermd.org/doc12/Amber15.pdf): If you have AmberTools15 and already `make install`, it is likely that you already have ``libsander``

# In[1]:


# load pytraj and load trajectories

import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

import pytraj as pt

traj = pt.iterload('tz2.nc', 'tz2.parm7')
print(traj)


# In[2]:


# compute the energy with igb=8 (GBneck2 solvation model)
data = pt.energy_decomposition(traj, igb=8)
data


# In[3]:


# get energy you want
print('potential energy', data['tot'])


# In[4]:


# solvation energy (igb=8)
print('solvation energy', data['gb'])


# In[5]:


# you can also get other energies.
# for example: data['dihedral']
list(data.keys())


# # Parallel calculation

# In[6]:


# serial: pt.energy_decomposition(traj, igb=8)

# parallel
data2 = pt.pmap(func=pt.energy_decomposition, traj=traj, igb=8, n_cores=4)
data2['gb']

