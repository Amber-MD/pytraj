from .parallel_mapping_mpi import pmap_mpi
from .parallel_mapping_multiprocessing import pmap
from pytraj.tools import concat_dict
from .pjob import PJob
from functools import partial
from pytraj import Frame


def _worker_state(rank, n_cores=1, traj=None, lines=[], dtype='dict', ref=None):
    # need to make a copy if lines since python's list is dangerous
    # it's easy to mess up with mutable list
    # do not use lines.copy() since this is not available in py2.7
    my_lines = [line for line in lines]
    from pytraj.utils import split_range
    from pytraj.core.cpp_core import _load_batch

    if ref is not None:
        if isinstance(ref, Frame):
            reflist = [ref, ]
        else:
            # list/tuplex
            reflist = ref
    else:
        reflist = []

    mylist = split_range(n_cores, 0, traj.n_frames)[rank]
    start, stop = mylist
    crdframes_string = 'crdframes ' + ','.join((str(start+1), str(stop)))

    for idx, line in enumerate(my_lines):
        if not line.lstrip().startswith('reference'):
            my_lines[idx] = ' '.join(('crdaction traj', line, crdframes_string))

    my_lines = ['loadtraj name traj',] + my_lines

    state = _load_batch(my_lines, traj)

    if reflist:
        for ref_ in reflist:
            ref_dset = state.data.add_new('reference')
            ref_dset.top = traj.top
            ref_dset.add_frame(ref_)
    state.run()
    for dset in state.data:
        if hasattr(dset, 'dtype') and dset.dtype == 'ref_frame':
            state.data.remove_set(dset)
    if dtype == 'dict':
        # exclude DatasetTopology and TrajectoryCpptraj
        return (rank, state.data[2:].to_dict())
    else:
        raise ValueError('must use dtype="dict"')


def _load_batch_pmap(n_cores=4, lines=[], traj=None, dtype='dict', root=0,
        mode='multiprocessing', ref=None):
    '''mpi or multiprocessing
    '''
    if mode == 'multiprocessing':
        from multiprocessing import Pool
        pfuncs = partial(_worker_state, n_cores=n_cores, traj=traj, dtype=dtype,
                lines=lines, ref=ref)
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
