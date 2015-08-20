
# coding: utf-8

# # **Aim**: download pdb file, calculate phi/psi for specific residue, show ramachandran plot

# In[1]:

get_ipython().magic('matplotlib inline')
from matplotlib import pyplot as plt
import seaborn as snb
snb.set()

import pytraj as pt


# In[2]:

# download trp-cage mini protein
# http://www.rcsb.org/pdb/explore.do?structureId=1l2y

traj = pt.load_pdb_rcsb('1L2Y')
print(traj)
print(set(res.name for res in traj.top.residues))


# In[3]:

# calculate phi/psi for Gly residues
# need to get indcies of Gly residues
indices = [idx for idx, res in enumerate(traj.top.residues) if 'GLY' in res.name]
print('Gly resdiue indices = ', indices)

ds = pt.multidihedral(traj, 'phi psi', resrange=indices)
print(ds)


# In[4]:

phi = pt.tools.flatten(ds.grep('phi'))
psi = pt.tools.flatten(ds.grep('psi'))

plt.xlim([-180, 180])
plt.ylim([-180, 180])
plt.xlabel('phi')
plt.ylabel('psi')
plt.plot(phi, psi, 'ro')

