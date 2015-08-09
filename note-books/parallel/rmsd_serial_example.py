"""calculat RMSD for 8 replica trajs using 1 cores.
Reference frame is the 1st frame of remd.x.000

System: 17443 atoms, 1000 frames, netcdf, 8 replicas (000 to 007), 200Mb/replica

python test_serial.py
"""

import numpy as np
from pytraj import io
import pytraj.common_actions as pyca

root_dir = "../../tests/data/nogit/remd/"
top_name = root_dir + "myparm.top"

# 8 replicas, 1000 frames
size = 8
sarr = np.empty((size, 1000))
REF = None
for i in range(size):
    fname = root_dir + "/remd.x.00" + str(i) # 000, 001, 002, 003 ...
    straj = io.load(fname, root_dir + "/myparm.parm7")
    if i == 0:
        REF = straj[0]
    sarr[i] = straj.calc_rmsd(REF, "@CA")
np.savetxt("./serial_rmsd.txt", sarr.flatten())
