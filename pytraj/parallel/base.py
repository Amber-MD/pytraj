import numpy as np
from functools import partial
from pytraj import Frame
from pytraj import pipe
from pytraj.datasets import CpptrajDatasetList
from pytraj.externals.six import string_types


def check_valid_command(commands):
    '''

    Parameters
    ----------
    commands : list/tuple of str

    Returns
    -------
    cmlist : newly updated command list
    need_ref : bool, whether to provide reference or not
    '''
    from pytraj.c_commands import analysis_commands

    if isinstance(commands, string_types):
        commands = [line.strip() for line in commands.split('\n') if line]
    else:
        commands = commands

    new_commands = commands[:]
    need_ref = False
    for idx, cm in enumerate(commands):
        cm = cm.strip()

        # try to insert 'refindex 0' if action requires ref but user does not provide
        if 'refindex' in cm:
            need_ref = True

        if ((cm.startswith('rms') or
             cm.startswith('nastruct') or
             cm.startswith('center') or
             cm.startswith('distrmsd') or
             cm.startswith('nativecontacts') or
             cm.startswith('symmetricrmsd')) and 'refindex' not in cm):

            cm = cm + ' refindex 0 '
            need_ref = True

        if cm.startswith('matrix'):
            raise ValueError('Not support matrix')

        for word in analysis_commands:
            if cm.startswith(word):
                raise ValueError(
                    'Not support cpptraj analysis keyword for parallel '
                    'calculation. You can use pmap for cpptraj actions to speed up the '
                    'IO and then perform '
                    'analysis in serial')
        new_commands[idx] = cm

    return new_commands, need_ref


def worker_byfunc(rank,
                  n_cores=None,
                  func=None,
                  traj=None,
                  args=None,
                  kwd=None,
                  iter_options=None,
                  apply=None):
    '''worker for pytraj's functions
    '''
    # need to unpack args and kwd
    if iter_options is None:
        iter_options = {}
    mask = iter_options.get('mask')
    rmsfit = iter_options.get('rmsfit')
    autoimage = iter_options.get('autoimage', False)
    iter_func = apply
    frame_indices = kwd.pop('frame_indices', None)

    if frame_indices is None:
        my_iter = traj._split_iterators(n_cores,
                                        rank=rank,
                                        mask=mask,
                                        rmsfit=rmsfit,
                                        autoimage=autoimage)
    else:
        my_indices = np.array_split(frame_indices, n_cores)[rank]
        my_iter = traj.iterframe(frame_indices=my_indices,
                                 mask=mask,
                                 rmsfit=rmsfit,
                                 autoimage=autoimage)
    n_frames = my_iter.n_frames
    kwd_cp = {}
    kwd_cp.update(kwd)

    if iter_func is not None:
        final_iter = iter_func(my_iter)
        kwd_cp['top'] = my_iter.top
    else:
        final_iter = my_iter

    data = func(final_iter, *args, **kwd_cp)
    return (rank, data, n_frames)


def worker_by_actlist(rank,
                      n_cores=2,
                      traj=None,
                      lines=None,
                      dtype='dict',
                      ref=None,
                      kwd=None):
    '''worker for cpptraj commands (string)
    '''
    # need to make a copy if lines since python's list is dangerous
    # it's easy to mess up with mutable list
    # do not use lines.copy() since this is not available in py2.7
    # Note: dtype is a dummy argument, it is always 'dict'
    if lines is None:
        lines = []
    frame_indices = kwd.pop('frame_indices', None)

    new_lines, need_ref = check_valid_command(lines)

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
        reflist = [traj[0],] if need_ref else []

    dslist = CpptrajDatasetList()

    if reflist:
        for ref_ in reflist:
            ref_dset = dslist.add('reference')
            ref_dset.top = traj.top
            ref_dset.add_frame(ref_)

    # create Frame generator
    fi = pipe(my_iter, commands=new_lines, dslist=dslist)

    # just iterate Frame to trigger calculation.
    for _ in fi:
        pass

    # remove ref
    return (rank, dslist[len(reflist):].to_dict())


def _load_batch_pmap(n_cores=4,
                     lines=None,
                     traj=None,
                     dtype='dict',
                     root=0,
                     mode='multiprocessing',
                     ref=None,
                     **kwd):
    '''mpi or multiprocessing
    '''
    if lines is None:
        lines = []
    if mode == 'multiprocessing':
        from multiprocessing import Pool
        pfuncs = partial(worker_by_actlist,
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
        data_chunk = worker_by_actlist(rank=rank,
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


def worker_state(rank, n_cores=1, traj=None, lines=None, dtype='dict'):
    '''worker for CpptrajState
    '''
    # need to make a copy if lines since python's list is dangerous
    # it's easy to mess up with mutable list
    # do not use lines.copy() since this is not available in py2.7
    if lines is None:
        lines = []
    my_lines = [line for line in lines]
    from pytraj.utils import split_range
    from pytraj.core.c_core import _load_batch

    mylist = split_range(n_cores, 0, traj.n_frames)[rank]
    start, stop = mylist
    crdframes_string = 'crdframes ' + ','.join((str(start + 1), str(stop)))

    for idx, line in enumerate(my_lines):
        if not line.lstrip().startswith('reference'):
            my_lines[idx] = ' '.join(('crdaction traj', line, crdframes_string
                                      ))

    # do not use 'extend' in this case
    # will get weird output
    my_lines = ['loadtraj name traj', ] + my_lines

    state = _load_batch(my_lines, traj)

    state.run()
    if dtype == 'dict':
        # exclude DatasetTopology and TrajectoryCpptraj
        return (rank, state.data[2:].to_dict())
    else:
        raise ValueError('must use dtype="dict"')
