import sys
import numpy as np
import pytraj as pt
import mdtraj as md
from pytraj.testing import Timer

fname, tname = ('GAAC3.1000frames.nc', 'GAAC.parm7')

def rmsd_():
    @Timer()
    def mdtraj_():
        print('mdtraj')
        traj = md.load(fname, top=tname)
        indices = traj.top.select('not water')
        return md.rmsd(traj, traj, 0, indices)
    
    @Timer()
    def pmap_(n_cores):
        print('n_cores = ', n_cores)
        traj = pt.iterload(fname, tname)
        print('traj size = %s (GB)' % traj._estimated_GB)
        x = pt.pmap(n_cores, pt.rmsd, traj, mask='!:WAT', ref=traj[0])
        return x

    print(mdtraj_())
    print(pmap_(n_cores=N_CORES))

if __name__ == '__main__':
    N_CORES = int(sys.argv[1])
    rmsd_()
