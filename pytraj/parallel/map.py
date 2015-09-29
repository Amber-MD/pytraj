import numpy as np
from pytraj.utils import split_range

def map(comm, calc_method, traj_or_list, command,
        root=0,
        dtype='ndarray', *args, **kwd):
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
        start, stop = split_range(size, 0, traj_or_list.n_frames)[rank]
        fa_chunk = traj_or_list(start=start, stop=stop) 
    dslist = calc_method(fa_chunk, command, dtype=dtype, *args, **kwd)
    total = comm.gather(dslist, root=root)
    return total
