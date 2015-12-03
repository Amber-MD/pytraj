import numpy as np
from functools import partial
from pytraj import Frame
from pytraj import create_pipeline
from pytraj.datasets import CpptrajDatasetList
from pytraj.externals.six import string_types


def check_valid_command(commands):
    '''

    Parameters
    ----------
    commands : list/tuple of str
    '''
    from pytraj.cpptraj_commands import analysis_commands

    if isinstance(commands, string_types):
        commands = [line.strip() for line in commands.split('\n') if line]
    else:
        commands = commands

    for cm in commands:
        cm = cm.strip()
        if cm.startswith('rms') and 'refindex' not in cm and 'reference' not in cm:
            raise ValueError('must prodive refindex for rms/rmsd command')
        if cm.startswith('matrix'):
            raise ValueError('Not support matrix')
        for word in analysis_commands:
            if cm.startswith(word):
                raise ValueError(
                    'Not support cpptraj analysis keyword for parallel '
                    'calculation. You can use pmap for cpptraj actions to speed up the '
                    'IO and then perform '
                    'analysis in serial')


def _worker_actlist(rank,
                    n_cores=2,
                    traj=None,
                    lines=[],
                    dtype='dict',
                    ref=None,
                    kwd=None):
    # need to make a copy if lines since python's list is dangerous
    # it's easy to mess up with mutable list
    # do not use lines.copy() since this is not available in py2.7
    frame_indices = kwd.pop('frame_indices', None)

    if frame_indices is None:
        my_iter = traj._split_iterators(n_cores, rank=rank)
    else:
        my_iter = traj.iterframe(
            frame_indices=np.array_split(frame_indices, n_cores)[rank])

    if ref is not None:
        if isinstance(ref, Frame):
            reflist = [ref, ]
        else:
            # list/tuplex
            reflist = ref
    else:
        reflist = []

    dslist = CpptrajDatasetList()

    if reflist:
        for ref_ in reflist:
            ref_dset = dslist.add_new('reference')
            ref_dset.top = traj.top
            ref_dset.add_frame(ref_)

    # create Frame generator
    fi = create_pipeline(my_iter, commands=lines, dslist=dslist)

    # just iterate Frame to trigger calculation.
    for _ in fi:
        pass

    # remove ref
    if dtype == 'dict':
        return (rank, dslist[len(reflist):].to_dict())
    else:
        raise ValueError('must use dtype="dict"')


def _load_batch_pmap(n_cores=4,
                     lines=[],
                     traj=None,
                     dtype='dict',
                     root=0,
                     mode='multiprocessing',
                     ref=None,
                     **kwd):
    '''mpi or multiprocessing
    '''
    if mode == 'multiprocessing':
        from multiprocessing import Pool
        pfuncs = partial(_worker_actlist,
                         n_cores=n_cores,
                         traj=traj,
                         dtype=dtype,
                         lines=lines,
                         ref=ref,
                         kwd=kwd)
        pool = Pool(n_cores)
        data = pool.map(pfuncs, range(n_cores))
        pool.close()
        pool.join()
        return data
    elif mode == 'mpi':
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.rank
        data_chunk = _worker_actlist(rank=rank,
                                     n_cores=n_cores,
                                     traj=traj,
                                     dtype=dtype,
                                     lines=lines,
                                     ref=ref,
                                     kwd=kwd)
        # it's ok to use python level `gather` method since we only do this once
        # only gather data to root, other cores get None
        data = comm.gather(data_chunk, root=root)
        return data
    else:
        raise ValueError('only support multiprocessing or mpi')
