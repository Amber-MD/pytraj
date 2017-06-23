
# coding: utf-8

# ## Using multiprocessing
# 
# ** Task **
# ```
# Perform RMSD calculation in parallel
# ```

# In[1]:

import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)
import pytraj as pt

traj = pt.iterload('tz2.nc', 'tz2.parm7')
print(traj)


# In[22]:

pt.pmap(pt.rmsd, traj, mask='@CA', ref=traj[0], n_cores=4)


# ## Using MPI
# 
# ** Task **
# ```
# Perform RMSD calculation in parallel
# ```

# In[24]:

get_ipython().run_cell_magic(u'file', u'my_script.py', u"\n## create a file name my_scrip.py (you can use any name you like)\n\nimport pytraj as pt\nfrom mpi4py import MPI\n\ncomm = MPI.COMM_WORLD\nrank = comm.rank\n\n# load files\ntraj = pt.iterload('tz2.nc', 'tz2.parm7')\n\n# call pmap_mpi for MPI\n\n# we dont need to specify n_cores=6 here since we will use `mpirun -n 6`\ndata = pt.pmap_mpi(pt.rmsd, traj, mask='@CA', ref=traj[0])\n\n# data is sent to rank==0\nif rank == 0:\n    print(data)\n\n# run")


# In[27]:

# run in shell

get_ipython().system(u' mpirun -n 4 python my_script.py')

