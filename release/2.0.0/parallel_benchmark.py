
# coding: utf-8

# ## System
# 
# ```python
# pytraj.TrajectoryIterator, 200000 frames: 
# Size: 58.150291 (GB)
# <Topology: 13008 atoms, 4189 residues, 4174 mols, PBC with box type = truncoct>
# 
# format: netcdf
# ```

# ## Multiprocessing

# <img src="./images/pmap_scaling/bench_pmap_casegroup.png">

# ## Distributed: using MPI (via mpi4py)

# <img src="./images/pmap_scaling/bench_pmap_mpi.png">

# ### script for multiprocessing
# 
# ```python
# from multiprocessing import cpu_count
# from glob import glob
# import pytraj as pt
# from time import time
# 
# print("max cores = {}".format(cpu_count()))
# filenames = glob('mdx/md*.nc') * 10
# traj = pt.iterload(filenames, 'prmtop')
# print(traj)
# 
# mask = '!:WAT,Na+'
# 
# t0 = time()
# pt.rmsd(traj, mask=mask)
# t_serial = time() - t0
# # print('serial time = ', t_serial)
# 
# func = pt.rmsd
# 
# # for n_cores in [1, 2, 4, 6, 8, 16, 20, 24]:
# for n_cores in [1, 2, 4, 6, 7]:
#     t0 = time()
#     pt.pmap(func, traj, mask=mask, n_cores=n_cores, ref=traj[0])
#     t_par = time() - t0
#     print(n_cores, t_serial / t_par)
# ```
# 

# ### script for MPI
# 
# ```python
# from glob import glob
# import pytraj as pt
# from time import time
# from mpi4py import MPI
# 
# comm = MPI.COMM_WORLD
# 
# # need to run program in serial and update serial time
# serial_time = 45.
# 
# filenames = glob('mdx/md*.nc') * 10
# traj = pt.iterload(filenames, 'prmtop')
# 
# mask = '!:WAT,Na+'
# 
# func = pt.rmsd
# 
# t0 = time()
# x = pt.pmap_mpi(func, traj, mask=mask, ref=traj[0])
# t_par = time() - t0
# 
# if comm.rank == 0:
#     print(serial_time/t_par)
# ```
