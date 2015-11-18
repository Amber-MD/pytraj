
# coding: utf-8

# In[1]:

# it's simple as long as you have internet and you're allowed to download

import pytraj as pt

traj = pt.io.loadpdb_rcsb("2KOC")
print (traj)

print ("coords of firs 10 atoms of 1st frame from pytraj")
print (traj[0, :10])

print ()
print ("calc_rmsd")
# print (traj.calc_rmsd(ref=0)) # all atoms, ref=first_frame
# use numpy to have prettier print
import numpy as np
import pytraj.common_actions as pyca
arr0 = np.asarray(pyca.rmsd(traj, ref=traj[0]))
print (arr0)


# # coords of firs 10 atoms of 1st frame from rbsb
# # http://www.rcsb.org/pdb/files/2KOC.pdb
# ```
# ATOM      1  OP3   G A   1      -8.886  -5.163   9.647  1.00 10.00           O  
# ATOM      2  P     G A   1     -10.305  -4.745  10.229  1.00 10.00           P  
# ATOM      3  OP1   G A   1     -10.282  -3.296  10.528  1.00 10.00           O  
# ATOM      4  OP2   G A   1     -10.676  -5.703  11.297  1.00 10.00           O  
# ATOM      5  O5'   G A   1     -11.278  -4.981   8.995  1.00 10.00           O  
# ATOM      6  C5'   G A   1     -11.721  -6.288   8.650  1.00 10.00           C  
# ATOM      7  C4'   G A   1     -12.147  -6.326   7.203  1.00 10.00           C  
# ATOM      8  O4'   G A   1     -12.908  -5.127   6.888  1.00 10.00           O  
# ATOM      9  C3'   G A   1     -11.018  -6.335   6.182  1.00 10.00           C  
# ATOM     10  O3'   G A   1     -10.515  -7.653   5.962  1.00 10.00           O
# ```

# <img src='http://www.rcsb.org/pdb/images/2koc_asym_r_500.jpg'>
