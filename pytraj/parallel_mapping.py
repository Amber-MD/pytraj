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
from pytraj import ired_vector_and_matrix, rotation_matrix
from pytraj import NH_order_parameters
from multiprocessing import cpu_count


def _concat_dict(iterables):
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


def _worker(rank,
           n_cores=None,
           func=None,
           traj=None,
           args=None,
           kwd=None,
           iter_options={}):
    # need to unpack args and kwd
    mask = iter_options.get('mask', None)
    rmsfit = iter_options.get('rmsfit', None)
    autoimage = iter_options.get('autoimage', False)
    my_iter = traj._split_iterators(n_cores, rank=rank, mask=mask, rmsfit=rmsfit,
            autoimage=autoimage)
    data = func(my_iter, *args, **kwd)
    return (rank, data, my_iter.n_frames)


def _pmap(func, traj, *args, **kwd):
    '''use python's multiprocessing to accelerate calculation. Limited calculations.

    Parameters
    ----------
    func : a pytraj's methods or a list of string or simply as a cpptraj' text
    traj : pytraj.TrajectoryIterator
    n_cores : int, number of cores to be used, default 2. Specify n_cores=-1 to use all available cores
    iter_options : dict, default {}
        Specify trajectory iterating option. This will be done before calling ``func``.
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

    vs::

        In [1]: pt.pmap(4, pt.radgyr, traj, dtype='dict')
        Out[1]:
        OrderedDict([('RoG_00000',
                      array([ 18.91114428,  18.93654996,  18.84969884,  18.90449256,
                              18.8568644 ,  18.88917208,  18.9430491 ,  18.88878079,
                              18.91669565,  18.87069722]))])

    This is experimental method, you should expect its syntax, default output will be changed.
    
    Examples
    --------
    >>> import numpy as np
    >>> import pytraj as pt
    >>> traj = pt.load_sample_data('tz2')
    >>> data = pt.pmap(pt.radgyr, traj, n_cores=4)
    >>> data # doctest: +SKIP
    [(0, array([ 18.91114428,  18.93654996]), 2),
     (1, array([ 18.84969884,  18.90449256]), 2),
     (2, array([ 18.8568644 ,  18.88917208]), 2),
     (3, array([ 18.9430491 ,  18.88878079,  18.91669565,  18.87069722]), 4)]
    >>> # in most cases, you can follow below command to join the data
    >>> pt.tools.flatten([x[1] for x in data]) # doctest: +SKIP
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

    >>> # cpptraj command style
    >>> data = pt.pmap(['distance :3 :7', 'vector mask :3 :12'], traj, n_cores=4)

    >>> # use iter_options
    >>> iter_options = {'autoimage': True, 'rmsfit': (0, '@CA')}
    >>> data = pt.pmap(pt.mean_structure, traj, iter_options=iter_options) 

    See also
    --------
    pytraj.parallel.map_mpi
    '''
    from multiprocessing import Pool
    from pytraj import TrajectoryIterator

    if 'n_cores' in kwd.keys():
        n_cores = kwd['n_cores']
        kwd.pop('n_cores')
    else:
        # 2 cores
        n_cores = 2
    if n_cores <= 0:
        # use all available cores
        n_cores = cpu_count()

    if 'iter_options' in kwd.keys():
        iter_options = kwd['iter_options']
        kwd.pop('iter_options')
    else:
        iter_options = {}

    if isinstance(func, (list, tuple)):
        # assume using _load_batch_pmap
        from pytraj.parallel import _load_batch_pmap
        data = _load_batch_pmap(n_cores=n_cores, traj=traj, lines=func, dtype='dict', root=0, mode='multiprocessing')
        return _concat_dict((x[1] for x in data))
    else:
        if not callable(func):
            raise ValueError('must callable argument')
        # pytraj's method
        if not hasattr(func, '_is_parallelizable'):
            raise ValueError("this method does not support parallel")
        elif not func._is_parallelizable:
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

        pfuncs = partial(_worker,
                         n_cores=n_cores,
                         func=func,
                         traj=traj,
                         args=args,
                         kwd=kwd,
                         iter_options=iter_options)

        data = p.map(pfuncs, [rank for rank in range(n_cores)])
        p.close()

        if func in [matrix.dist, matrix.idea]:
            mat  = np.sum((val[1] * val[2] for val in data)) / traj.n_frames
            return mat
        elif func in [ired_vector_and_matrix, ]:
            # data is a list of (rank, (vectors, matrix), n_frames)
            mat = np.sum((val[1][1] * val[2] for val in data)) / traj.n_frames
            #vecs = np.vstack((val[1][0] for val in data))
            vecs = np.column_stack(val[1][0] for val in data)
            return (vecs, mat)
        elif func in [rotation_matrix, ]:
            if 'with_rmsd' in kwd.keys() and kwd['with_rmsd']:
                # data is a list of (rank, (mat, rmsd), n_frames)
                mat = np.row_stack(val[1][0] for val in data)
                rmsd_ = np.hstack(val[1][1] for val in data)
                return mat, rmsd_
            else:
                mat = np.row_stack(val[1] for val in data)
                return mat
        elif func == mean_structure:
            xyz = np.sum((x[2] * x[1].xyz for x in data)) / traj.n_frames
            frame = Frame(xyz.shape[0])
            frame.xyz[:] = xyz
            return frame
        else:
            if dtype == 'dict':
                new_dict = _concat_dict((x[1] for x in data))
                return new_dict
            else:
                return data

def pmap(func=None, traj=None, *args, **kwd):
    if func != NH_order_parameters:
        return _pmap(func, traj, *args, **kwd)
    else:
        if 'n_cores' in kwd.keys():
            if kwd['n_cores'] == 1:
                # use default n_cores=2 instead of 1
                kwd['n_cores'] = 2
            if kwd['n_cores'] <= 0:
                kwd['n_cores'] = cpu_count()
        else:
            # use n_cores=2 for default value
            kwd['n_cores'] = 2
        return NH_order_parameters(traj, *args, **kwd)

pmap.__doc__ = _pmap.__doc__
