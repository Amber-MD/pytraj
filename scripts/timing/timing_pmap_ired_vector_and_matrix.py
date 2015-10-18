import sys
import numpy as np
import pytraj as pt
from pytraj.testing import Timer
from pytraj.testing import aa_eq

traj = pt.iterload('traj.nc', 'lysozyme.top', frame_slice=(0, 10000))

h = traj.top.select('@H')
n = h - 1
nh = list(zip(n, h))

def ired_():
    @Timer()
    def normal_(traj=traj, nh=nh):
        print('serial')
        return pt.ired_vector_and_matrix(traj, nh)
    
    @Timer()
    def pmap_(traj=traj, nh=nh, n_cores=4):
        print('n_cores = ', n_cores)
        print(traj)
        return pt.pmap(pt.ired_vector_and_matrix, traj, nh, n_cores=n_cores)

    x0 = normal_()
    x1 = pmap_(n_cores=N_CORES)

    #aa_eq(x0[0], x1[0])
    #aa_eq(x0[1], x1[1])

if __name__ == '__main__':
    N_CORES = int(sys.argv[1])
    ired_()

# output:
# not bad
#    serial
#    15.39 (s)

#    n_cores =  6
#    2.62 (s)
