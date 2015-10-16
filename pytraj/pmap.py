from functools import partial

def worker(rank,
           n_cores=None,
           func=None,
           traj=None,
           args=None,
           kwd=None):
    # need to unpack args and kwd
    my_iter = traj._split_iterators(n_cores, rank=rank)
    data = func(my_iter, *args, **kwd)
    return (rank, data, my_iter.n_frames)


def pmap(n_cores=2, func=None, traj=None, *args, **kwd):
    '''

    Returns
    -------
    out : list of (rank, data)

    Notes
    -----
    If calculation require a reference structure, users need to explicit provide reference
    as a Frame (not an integer number)
    
    Examples
    --------
    >>> import numpy as np
    >>> import pytraj as pt
    >>> traj = pt.load_sample_data('tz2')
    >>> data = pt.pmap(4, pt.radgyr, traj=traj)
    >>> data
    [(0, array([ 18.91114428,  18.93654996])),
     (1, array([ 18.84969884,  18.90449256])),
     (2, array([ 18.8568644 ,  18.88917208])),
     (3, array([ 18.9430491 ,  18.88878079,  18.91669565,  18.87069722]))]
    >>> # in most cases, you can follow below command to join the data
    >>> pt.tools.flatten([x[1] for x in data])
    [18.911144277821389,
     18.936549957265814,
     18.849698842157373,
     18.904492557176411,
     18.856864395949234,
     18.889172079501037,
     18.943049101357886,
     18.888780788130308,
     18.916695652897396,
     18.870697222142766]
    '''
    from multiprocessing import Pool
    from pytraj import TrajectoryIterator

    if not hasattr(func, '_is_parallelizable') or not func._is_parallelizable:
        raise ValueError("this method does not support parallel")

    if not isinstance(traj, TrajectoryIterator):
        raise ValueError('only support TrajectoryIterator')

    p = Pool(n_cores)
    pfuncs = partial(worker,
                     n_cores=n_cores,
                     func=func,
                     traj=traj,
                     args=args,
                     kwd=kwd)
    result = p.map(pfuncs, [rank for rank in range(n_cores)])
    p.close()
    return result
