from pytraj.utils import _import_numpy

has_np, np = _import_numpy()


def map(comm, calc_method, traj_or_list, command, root=0, dtype='ndarray', *args, **kwd):
    """

    # creat file
    cat > test_mpi.py <<EOF
    from mpi4py import MPI
    from pytraj.parallel import map as pymap
    from pytraj import io
    import pytraj.common_actions as pyca

    comm = MPI.COMM_WORLD
    fa = io.load(traj_name, parm_name)[:] 
    if rank == 0:
        result_arr = pymap(comm, fa, pyca.calc_molsurf, "@CA")
        print (result_arr)
    EOF

    # run in parallel
    mpirun -n 4 python ./test_mpi.py  
    """
    size = comm.size
    rank = comm.rank
    # create iterator
    if isinstance(traj_or_list, (list, tuple)):
        # split traj_or_list into n_cores
        fa_chunk = traj_or_list[rank]
    else:
        # split single traj into n_cores
        traj = traj_or_list
        chunk = traj.n_frames // size

        if rank == size - 1:
            fa_chunk = traj(start=rank * chunk)
        else:
            fa_chunk = traj(start=rank * chunk, stop=(rank + 1) * chunk - 1)

    dslist = calc_method(fa_chunk, command, dtype=dtype, *args, **kwd)
    total = comm.gather(dslist, root=root)
    return total
