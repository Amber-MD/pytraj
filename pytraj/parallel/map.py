from pytraj.utils import _import_numpy

has_np, np = _import_numpy()

def map(comm, farray, calc_method, command, *args, **kwd):
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
    chunk = farray.n_frames // size
    fa_chunk = farray[rank*chunk : (rank + 1) *chunk]
    dslist = calc_method(command, fa_chunk, *args, **kwd)
    arr0 = dslist.to_ndarray()
    return arr0
