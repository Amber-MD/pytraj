"""calculat RMSD for 8 replica trajs using openmp with 8 cores
Reference frame is the 1st frame of remd.x.000

System: 17443 atoms, 1000 frames, netcdf, 8 replicas (000 to 007), 200Mb/replica

python test_openmp_mdtraj.py
"""

import numpy as np
import mdtraj as md

size = 8
sarr = np.empty((size, 1000))
REF = None

root_dir = "../../tests/data/nogit/remd/"

for i in range(size):
    fname = root_dir + "/remd.x.00" + str(i)
    straj = md.load_netcdf(fname, root_dir + "/myparm.parm7")
    indices = straj.top.select("name CA")
    if i == 0:
        REF = straj[0]
    sarr[i] = md.rmsd(straj, REF, 0, indices)
np.savetxt("rmsd_mdtraj_openmp.txt", sarr.flatten())
