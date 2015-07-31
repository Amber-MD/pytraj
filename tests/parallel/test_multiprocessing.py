import pytraj as pt
import sys

trajname = './data/nogit/remd/remd.x.000'
topname = './data/nogit/remd/myparm.top'
traj = pt.iterload(trajname, topname)

n_jobs = int(sys.argv[1])

def worker(i):
    traj = pt.iterload(trajname, topname)
    print(pt.radgyr(traj(stop=1000)))

if sys.argv[2] == 'p': 
    print('parallel')
    from multiprocessing import Pool

    p = Pool(n_jobs)
    p.map(worker, [_ for _ in range(n_jobs)])

elif sys.argv[2] == 's': 
    print('serial')
    for i in range(n_jobs):
        worker(i)

else:
    raise ValueError('must have sys.argv[2]')
