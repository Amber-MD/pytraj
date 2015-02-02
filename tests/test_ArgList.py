from pytraj.base import *
from pytraj.TrajectoryIO import TrajectoryIO

arglist = """
parm Tc5b.top
reference Tc5b.nat.crd
reftraj ./md1_prod.x parmindex 100
trajin ../../output_md1_prod/md1_prod.x 1 20000
dihedral psi1 :1@N :1@CA :1@C :2@N out test1.dat
dihedral phi2 :1@C :2@N :2@CA :2@C out PHIPSI_time_psi1_phi2_psi2_phi3.dat
dihedral psi2 :2@N :2@CA :2@C :3@N out test2.dat
dihedral phi3 :2@C :3@N :3@CA :3@C out test3.dat
rmsd reftraj ./md1_prod.x parmindex 100 out rmsd.dat :3-18@CA
analyze 
go"""

actionArgs = ArgList(arglist)
print(actionArgs.get_string_key("out"))
print(actionArgs.get_string_key("out"))
print(actionArgs.get_string_key("out"))
print(actionArgs.get_string_key("out"))

traj = Trajectory()
print(traj.check_frame_args(ArgList("1 100 10"), 9999))
print(traj.setup_frame_info()) 

top = Topology("data/Tc5b.top")
#trajio = TrajectoryIO()
#print trajio.total_frames("./data/md1_prod.Tc5b.x", top)
