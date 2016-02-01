import numpy as np
from pytraj.utils import split_range
from pytraj.tools import concat_dict
from pytraj.get_common_objects import get_reference


def pmap_mpi(func, traj, *args, **kwd):
    """parallel with MPI (mpi4py)

    Parameters
    ----------
    func : a function
    traj : pytraj.TrajectoryIterator
    *args, **kwd: additional arguments

    Examples
    --------
    .. code-block:: bash

        $ # create test_radgyr.py file
        $ cat > test_radgyr.py <<EOF
        import pytraj as pt
        from mpi4py import MPI
        comm = MPI.COMM_WORLD

        traj = pt.iterload('tz2.nc', 'tz2.parm7')

        result_arr = pt.pmap_mpi(pt.radgyr, traj, "@CA")

        if comm.rank == 0:
            # save data to disk to read later by pytraj.read_pickle
            # pt.to_pickle(result_arr, 'output.pk')
            print(result_arr)
        EOF

        $ # run in parallel
        $ mpirun -n 4 python ./test_radgyr.py
        [array([ 8.10916061,  7.7643485 ,  8.09693108, ...,  9.70825678,
                9.3161563 ,  8.86720964]), array([ 8.82037273,  8.89008289,  9.48540176, ...,  9.29585981,
                9.53138062,  9.19155977]), array([ 9.13735723,  8.94651001,  8.97810478, ...,  7.68751186,
                8.31361647,  7.83763754]), array([ 7.37423766,  7.05637263,  6.52135566, ...,  6.38061648,
                6.24139008,  6.48994552])]
    """
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    n_cores = comm.size
    rank = comm.rank

    # update reference
    if 'ref' in kwd:
        kwd['ref'] = get_reference(traj, kwd['ref'])

    if not isinstance(func, (list, tuple)):
        # split traj to ``n_cores`` chunks, perform calculation
        # for rank-th chunk
        if 'dtype' not in kwd:
            kwd['dtype'] = 'dict'

        frame_indices = kwd.pop('frame_indices', None)
        if frame_indices is None:
            start, stop = split_range(n_cores, 0, traj.n_frames)[rank]
            my_iter = traj.iterframe(start=start, stop=stop)
        else:
            my_indices = np.array_split(frame_indices, n_cores)[rank]
            my_iter = traj.iterframe(frame_indices=my_indices)
        data = func(my_iter, *args, **kwd)
        total = comm.gather(data, root=0)
        if rank == 0:
            total = concat_dict(x for x in total)
    else:
        # cpptraj command style
        from pytraj.parallel.base import _load_batch_pmap
        total = _load_batch_pmap(n_cores=n_cores,
                                 traj=traj,
                                 lines=func,
                                 dtype='dict',
                                 root=0,
                                 mode='mpi',
                                 **kwd)
        if rank == 0:
            # otherwise, total=None
            total = concat_dict((x[1] for x in total))
    return total
