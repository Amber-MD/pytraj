import numpy as np
from pytraj.utils import split_range

def map_mpi(func, traj, command, dtype='ndarray', *args, **kwd):
    """parallel with MPI (mpi4py)

    Paramters
    ---------
    func : a function
    traj : pytraj.TrajectoryIterator
    command : str
    dtype : default 'ndarray'
        dtype of return output
    *args, **kwd: additional arguments

    Examples
    --------
    # creat ``test_mpi.py`` file
    cat > test_mpi.py <<EOF
    import pytraj as pt
    from pytraj import io
    from pytraj.parallel import map as pymap

    fa = pt.iterload(traj_name, parm_name)
    if rank == 0:
        result_arr = pymap(pyca.calc_molsurf, traj, "@CA")
        print (result_arr)
    EOF

    # run in parallel
    mpirun -n 4 python ./test_mpi.py  
    """
    from mpi4py import MPI
    comm = MPI.COMM_WORLD 
    size = comm.size
    rank = comm.rank

    # split traj to ``size`` chunks, perform calculation 
    # for rank-th chunk
    start, stop = split_range(size, 0, traj.n_frames)[rank]
    fa_chunk = traj(start=start, stop=stop) 

    dslist = func(fa_chunk, command, dtype=dtype, *args, **kwd)

    # gather data to root
    total = comm.gather(dslist, root=0)
    return total

map = map_mpi 
