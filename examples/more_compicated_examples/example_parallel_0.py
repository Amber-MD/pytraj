# parallel: 4 cores
# vs serial vs mdtraj
import sys
from glob import glob
import pytraj as pt
try:
    import mdtraj as md
except ImportError:
    print("no mdtraj. skip")

flist = sorted(glob('../data/nogit/remd/remd.x.*'))[:8]
print(flist)
top = pt.load_topology('../data/nogit/remd/myparm.parm7')

if sys.argv[1] == 'parallel':
    print('parallel')
    from multiprocessing import Pool
    top = pt.load_topology('../data/nogit/remd/myparm.parm7')

    def worker(rank, n_cores=4):
        traj = pt.iterload(flist, top)
        return (rank, pt.calc_radgyr(
            traj.split_iterators(n_cores,
                                 rank=rank),
            top=traj.top))

    p = Pool(4)
    result = p.map(worker, [rank for rank in range(4)])
    print(result)

elif sys.argv[1] == 'serial':
    top = pt.load_topology('../data/nogit/remd/myparm.parm7')
    print('serial')
    traj = pt.iterload(flist, top=top)
    print(pt.calc_radgyr(traj))

elif sys.argv[1] == 'mdtraj':
    print('mdtraj')
    m_top = md.load_topology(top.filename)
    for fname in flist:
        m_traj = md.load_netcdf(fname, top=m_top)
        print(md.compute_rg(m_traj))

else:
    raise ValueError('must have sys.argv[2]')
