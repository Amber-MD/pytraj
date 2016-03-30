import sys
import numpy as np
import pytraj as pt
from pytraj.testing import Timer

traj = pt.iterload('GAAC3.5000frames.nc', 'GAAC.topo', frame_slice=(0, 5000))


def distance_molsurf():
    # distance, molsurf
    @Timer()
    def normal_(traj=traj):
        print('serial')
        state = pt.load_batch(traj, '''
        distance :3 :2
        molsurf !:WAT
        ''')
        state.run()
        return state.data[1:].values
    
    @Timer()
    def pmap_(traj=traj, n_cores=4):
        print('n_cores = ', n_cores)
        x = pt._load_batch_pmap(n_cores, traj, lines=['distance :3 :2', 'molsurf !:WAT'])
        return x
    
    normal_()
    pmap_(n_cores=N_CORES)

def matrix_():
    @Timer()
    def normal_(traj=traj):
        print('serial')
        state = pt.load_batch(traj, '''
        matrix dist @P
        ''')
        state.run()
        return state.data[1:].values
    
    @Timer()
    def pmap_(traj=traj, n_cores=4):
        print('n_cores = ', n_cores)
        print(traj)
        x = pt.pmap(n_cores, pt.matrix.dist, traj, '@P', dtype='ndarray')
        # make sure to reproduce serial version
        #data = np.sum((v[1] * v[2] for  v in x), axis=0)
        #return data / traj.n_frames
        return x

    #print(normal_())
    print(pmap_(n_cores=N_CORES))

def rmsd_():
    @Timer()
    def normal_(traj=traj):
        print('serial')
        state = pt.load_batch(traj, '''
        rms !:WAT
        ''')
        state.run()
        return state.data[1:].values
    
    @Timer()
    def pmap_(traj=traj, n_cores=4):
        print('n_cores = ', n_cores)
        print(traj)
        x = pt.pmap(n_cores, pt.rmsd, traj, mask='!:WAT', ref=traj[0])
        return x

    print(normal_())
    print(pmap_(n_cores=N_CORES))

if __name__ == '__main__':
    N_CORES = int(sys.argv[1])
    #matrix_()
    rmsd_()
