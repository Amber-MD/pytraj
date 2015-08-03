from functools import partial

def worker(rank, n_cores=None, func=None, traj=None, *args, **kwd):
    from pytraj import iterload
    local_traj = iterload(traj.filelist, top=traj.top)
    return (rank, func(
            traj.split_iterators(n_cores, rank=rank),
            top=traj.top, *args, **kwd))

def pmap(n_cores=2, func=None, traj=None, *args, **kwd):
    print(traj)
    from multiprocessing import Pool

    p = Pool(n_cores)
    pfuncs = partial(worker, n_cores=n_cores,
                     func=func, traj=traj, args=args, kwd=kwd)
    return p.map(pfuncs, [rank for rank in range(n_cores)])
