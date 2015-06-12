from pytraj import Trajectory, io

traj = Trajectory(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
traj2 = io.iterload('./data/DPDP.nc', "./data/DPDP.parm7")
top = traj.top
frame = traj[0]
