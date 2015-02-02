import numpy as np
from pycpptraj import *
from pycpptraj.io import writeparm, writetraj


act = allactions.Action_RandomizeIons()

traj = load(filename="../adh206.tip3p.rst7.gz", top="../adh206.ff10.tip3p.parm7.gz")
print traj.top.n_atoms
print traj.size
farray = FrameArray()

for i in range(10):
    farray.append(traj[0])

farray.top = traj.top.copy()
newframe = Frame()
oldframe = farray[0].copy()
print oldframe.coords == farray[0].coords

act.master(command="randomizeions @Na+ around :1-16 by 5.0 overlap 3.0", 
           currenttop=traj.top, 
           currentframe=farray[0],
           )

print farray.top.n_atoms
print oldframe.coords == farray[0].coords

#writetraj(filename="out.crd", traj=farray[0], top=traj.top, fmt='amberrestart', overwrite=True)
newtraj = load(filename="../random.crd.save", top=traj.top)
savedframe = newtraj[0]

print savedframe.rmsd(oldframe)
print savedframe.rmsd(farray[0])
