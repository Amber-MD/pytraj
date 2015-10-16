# do not use relative import here. Treat this module as a seperated package.
from functools import partial
from pytraj.cpp_options import info as compiled_info
from collections import OrderedDict
import numpy as np
from pytraj.externals.six import string_types, iteritems
from pytraj.datasetlist import stack
from pytraj._get_common_objects import _get_data_from_dtype
from pytraj import matrix
from pytraj import mean_structure 
from pytraj import Frame


def concat_dict(iterables):
    # we have this function in pytraj.tools but copy here to be used as internal method
    # TODO: fill missing values?
    """concat dict

    iterables : iterables that produces OrderedDict
    """
    new_dict = OrderedDict()
    for i, d in enumerate(iterables):
        if i == 0:
            # make a copy of first dict
            new_dict.update(d)
        else:
            for k, v in iteritems(new_dict):
                new_dict[k] = np.concatenate((new_dict[k], d[k]))
    return new_dict


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
    '''use python's multiprocessing to accelerate calculation. Limited calculations.

    Parameters
    ----------
    n_cores : int, number of cores to be used, default 2
    func : a pytraj's methods or a list of string or simply as a cpptraj' text
    traj : pytraj.TrajectoryIterator
    *args, **kwd: additional keywords

    Returns
    -------
    out : if dtype='dict', return an OrderedDict, data is automatically joint. if dtype is
          not 'dict', return a list of (rank, data, n_frames), data is NOT automatically
          joint. Note that for [matrix.dist, matrix.idea, mean_structure] calculation,
          data is always joint (does not depend on dtype)

    Notes
    -----
    If calculation require a reference structure, users need to explicit provide reference
    as a Frame (not an integer number). For example, pt.pmap(4, pt.rmsd, traj, ref=-3)
    won't work, use ``ref=traj[3]`` instead.

    This method only benifits you if your calculation is quite long (saying few minutes to
    few hours). For calculation that takes less than 1 minutes, you won't see the
    significant speed up (or even slower) since pytraj need to warm up and need to gather
    data when the calculation done.

    The parallel cacluation is very simple, trajectory will be split (almost equal) to
    different chunk (n_chunks = n_cores), pytraj/cpptraj perform calculation for each
    chunk in each core and then send data back to master. Note that we are using Python's
    built-in multiprocessing module, so you can use this method interactively in Ipython
    and ipython/jupyter notebook. This behavior is different from using MPI, in which you
    need to write a script, escaping ipython ession and type something like::
        
        mpirun -n 4 python my_script.py
    
    Examples
    --------
    >>> import numpy as np
    >>> import pytraj as pt
    >>> traj = pt.load_sample_data('tz2')
    >>> data = pt.pmap(4, pt.radgyr, traj=traj)
    >>> data
    [(0, array([ 18.91114428,  18.93654996]), 2),
     (1, array([ 18.84969884,  18.90449256]), 2),
     (2, array([ 18.8568644 ,  18.88917208]), 2),
     (3, array([ 18.9430491 ,  18.88878079,  18.91669565,  18.87069722]), 4)]
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

    if isinstance(func, (list, tuple, string_types)):
        # assume using _load_batch_pmap
        from pytraj.parallel import _load_batch_pmap
        data = _load_batch_pmap(n_cores=n_cores, traj=traj, lines=func, dtype='dict', root=0, mode='multiprocessing')
        return concat_dict((x[1] for x in data))
    else:
        # pytraj's method
        if not hasattr(func, '_is_parallelizable') or not func._is_parallelizable:
            raise ValueError("this method does not support parallel")
        else:
            if hasattr(func, '_openmp_capability') and func._openmp_capability and 'OPENMP' in compiled_info():
                raise RuntimeError("this method supports both openmp and pmap, but your cpptraj "
                "version was installed with openpm. Should not use both openmp and pmap at the "
                "same time. In this case, do not use pmap since openmp is more efficient")

        if not isinstance(traj, TrajectoryIterator):
            raise ValueError('only support TrajectoryIterator')

        p = Pool(n_cores)
        if 'dtype' in kwd.keys():
            dtype = kwd['dtype']
        else:
            dtype = None

        pfuncs = partial(worker,
                         n_cores=n_cores,
                         func=func,
                         traj=traj,
                         args=args,
                         kwd=kwd)

        data = p.map(pfuncs, [rank for rank in range(n_cores)])
        p.close()

        if func in [matrix.dist, matrix.idea]:
            y = np.sum((val[1] * val[2] for val in data)) / traj.n_frames
            return y
        elif func == mean_structure:
            xyz = np.sum((x[2] * x[1].xyz for x in data)) / traj.n_frames
            frame = Frame(xyz.shape[0])
            frame.xyz[:] = xyz
            return frame
        else:
            if dtype == 'dict':
                new_dict = concat_dict((x[1] for x in data))
                return new_dict
            else:
                return data
