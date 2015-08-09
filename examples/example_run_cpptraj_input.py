import pytraj as pt
from pytraj.datafiles import load_cpptraj_output
import os

# use iterload for memory saving
traj = pt.iterload("../tests/data/md1_prod.Tc5b.x", "../tests/data/Tc5b.top")

trajin = """
parm ../tests/data/Tc5b.top
trajin ../tests/data/md1_prod.Tc5b.x
rms first @CA
distance :2 :3
matrix dist
"""

cout = load_cpptraj_output(trajin)
print(cout)
