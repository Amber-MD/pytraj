import pytraj as pt
from functools import partial
from glob import glob
import sys

top = pt.load_topology('../data/nogit/remd/myparm.top')
flist = sorted(glob('../data/nogit/remd/remd.x.*')[:40])

if sys.argv[1] == 'parallel': 
    print('parallel')
    from multiprocessing import Pool

    def worker(rank, n_cores=4):
        traj = pt.iterload(flist, top)
        print("MB = ", traj._estimated_MB)
        print(traj)
        print('n_frames = ', pt.calc_radgyr(
             traj.split_iterators(
                 n_cores, rank=rank), top=traj.top).shape)

    p = Pool(4)
    p.map(worker, [rank for rank in range(4)])

elif sys.argv[1] == 'serial': 
    print('serial')
    traj = pt.iterload(flist, top=top)
    print("MB = ", traj._estimated_MB)
    print(traj)
    print('n_frames = ', pt.calc_radgyr(traj).shape)

else:
    raise ValueError('must have sys.argv[2]')
