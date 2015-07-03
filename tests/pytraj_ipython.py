from pytraj import Trajectory, io
import pytraj as pt

traj = Trajectory(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
traj2 = io.iterload('./data/DPDP.nc', "./data/DPDP.parm7")
top = traj.top
frame = traj[0]

dslist = traj.calc_multidihedral()
