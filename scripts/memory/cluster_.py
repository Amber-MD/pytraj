import pytraj as pt
from memory_profiler import profile

mask = '!:WAT,K+,Cl-'

@profile
def calc_():
    traj = pt.iterload('GAAC3.5000frames.nc', 'GAAC.topo', frame_slice=(0, 1000))
    print(traj)
    pt.cluster.kmeans(traj, mask=mask, n_clusters=10)

calc_()
