# require: mpi4py
# Example: mpirun -n 4 mpi_pysander.py

# import MPI to get rank
from mpi4py import MPI
import pytraj as pt
from pytraj.testing import aa_eq

try:
    import sander
    comm = MPI.COMM_WORLD
    # end. you are free to update anything below here

    # split remd.x.000 to N cores and do calc_surf in parallel
    root_dir = "../../tests/data/"
    traj_name = root_dir + "tz2.nc"
    parm_name = root_dir + "tz2.parm7"

    # load to TrajectoryIterator
    traj = pt.iterload(traj_name, parm_name)
    inp = sander.gas_input(8)

    # gather the data
    # if rank != 0: data is None
    data = pt.pmap_mpi(pt.energy_decomposition, traj, mm_options=inp)

    if comm.rank == 0:
        # make sure to reproduce serial output
        serial_data = pt.energy_decomposition(traj, mm_options=inp)
        aa_eq(pt.tools.dict_to_ndarray(data), pt.tools.dict_to_ndarray(serial_data))
except ImportError:
    print('does not have sander. skip this example')
