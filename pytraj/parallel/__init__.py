from .map import map_mpi
from pytraj.tools import concat_dict
from .pjob import PJob
from functools import partial


def get_comm_size_rank():
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.rank
        size = comm.size
        return comm, size, rank
    except ImportError:
        comm = rank = size = None


def gather(name='data', clients=None, restype='ndarray'):
    '''gather data from different clients

    Parameters
    ----------
    name : name of the output holder
        for example: data = pytraj.calc_radgyr(traj) --> name = 'data'
    clients : IPython.parallel.Client objects
        number of clients == n_cores you use
    restype : str, {'ndarray', 'dataset'}, default 'ndarray' 
        if 'ndarray': hstack data by numpy.vstack
        if 'dataset': 'data' should be a list of dict, then will be converted
        to `pytraj.datasetlist.DatasetList` object

    Examples
    --------
    (fill me)
    '''
    if restype == 'ndarray':
        import numpy as np
        return np.hstack((x[name] for x in clients))
    elif restype == 'dataset':
        # it's user's responsibility to return a list of dicts
        from pytraj import datasetlist
        iter_of_dslist = (
            datasetlist._from_full_dict(x[name]) for x in clients)
        return datasetlist.vstack(iter_of_dslist)
    elif restype == 'dict':
        return concat_dict((x[name] for x in clients))
    else:
        raise ValueError("must be ndarray | dataset | dict")


def _worker_state(rank, n_cores=1, traj=None, lines=[], dtype='dict'):
    # need to make a copy if lines since python's list is dangerous
    # it's easy to mess up with mutable list
    # do not use lines.copy() since this is not available in py2.7
    my_lines = [line for line in lines]
    from pytraj.utils import split_range
    from pytraj.core.cpp_core import _load_batch

    mylist = split_range(n_cores, 0, traj.n_frames)[rank]
    start, stop = mylist
    crdframes_string = 'crdframes ' + ','.join((str(start+1), str(stop)))

    for idx, line in enumerate(my_lines):
        if not line.lstrip().startswith('reference'):
            my_lines[idx] = ' '.join(('crdaction traj', line, crdframes_string))

    my_lines = ['loadtraj name traj',] + my_lines

    state = _load_batch(my_lines, traj)
    state.run()
    if dtype == 'dict':
        # exclude DatasetTopology and TrajectoryCpptraj
        return (rank, state.data[2:].to_dict())
    elif dtype == 'state':
        return state

def _load_batch_pmap(n_cores=4, traj=None, lines=[], dtype='dict', root=0, mode='multiprocessing'):
    '''mpi or multiprocessing
    '''
    if mode == 'multiprocessing':
        from multiprocessing import Pool
        pfuncs = partial(_worker_state, n_cores=n_cores, traj=traj, dtype=dtype, lines=lines)
        pool = Pool(n_cores)
        data = pool.map(pfuncs, range(n_cores))
        pool.close()
        pool.join()
        return data
    elif mode == 'mpi':
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        size = comm.size
        rank = comm.rank
        data_chunk = _worker_state(rank, n_cores=size, traj=traj, lines=lines, dtype=dtype)
        # it's ok to use python level `gather` method since we only do this once
        # only gather data to root, other cores get None
        data = comm.gather(data_chunk, root=root)
        return data
    else:
        raise ValueError('only support multiprocessing or mpi')
