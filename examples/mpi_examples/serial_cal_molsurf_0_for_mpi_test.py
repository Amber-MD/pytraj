# (require: mpi4py, numpy)
# python serial_cal_molsurf_0.py

from pytraj import io
import pytraj.common_actions as pyca

root_dir = "../../tests/data/nogit/remd/"
traj_name = root_dir + "/remd.x.000"
parm_name = root_dir + "myparm.top"

traj = io.load(traj_name, parm_name)
arr = pyca.calc_molsurf(traj, "@CA")
print(arr[:10])
