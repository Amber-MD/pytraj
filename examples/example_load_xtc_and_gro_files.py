import pytraj as pt
import mdtraj as md

traj = pt.load("./nvt.xtc", "nvt.gro", engine='mdtraj')
