# datafiles and this test.py were provided by jmborr
# issue #991: # Error "ValueError: Buffer not C contiguous" 
# https://github.com/Amber-MD/pytraj/issues/991

import os
import argparse
import numpy
import pytraj
from pdb import set_trace as tr

traj = pytraj.load('short.dcd', 'pdb', mask='(!:1-256)&(@H1,@H2,@H3,@H4,@H5)')

# Trajectory of the center of mass of the first two residues
minitraj =  numpy.empty((2,traj.n_frames,3))
minitraj[0] = pytraj.center_of_mass(traj,mask=':1')
minitraj[1] = pytraj.center_of_mass(traj,mask=':2')
minitraj = minitraj.transpose((1,0,2))

# "topology" using the first two atoms
minitop = traj.top['@1,2']

# Save trajectory
pytraj.write_traj(filename='/tmp/junk.crd', traj=minitraj, top=minitop, overwrite=True)
