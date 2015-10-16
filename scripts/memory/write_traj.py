from pytraj.sandbox import write_traj
import pytraj as pt
from memory_profiler import profile

@profile
def test():
    traj = pt.iterload('GAAC3.5000frames.nc', 'GAAC.topo')
    print(traj._estimated_GB)
    write_traj('test.nc', traj, frame_indices=range(500))

test()
