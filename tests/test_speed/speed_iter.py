#!/usr/bin/env python
import pytraj as pt
from pytraj.testing import Timer, aa_eq
import numpy as np
import mdtraj as md
'''compare iteration speed between pytraj.Trajectory and cpptraj's CRD set
'''

n_frames = 100
traj = pt.iterload('md.nc', 'tc5bwat.parm7', frame_slice=(0, n_frames))
m_traj = md.load(traj.filename, top=traj.top.filename, indices=range(n_frames))
print(m_traj)


@Timer()
def iter_(traj):
    print(traj)
    for _ in traj:
        for __ in traj:
            _.rmsd(__)


@Timer()
def rms2d_():
    state = pt.load_cpptraj_state('''
            parm {0}
            loadcrd md.trj 0 {1} crd_
            #crdaction crd_ rms2d
            rms2d crdset crd_ out test.dat
            '''.format(traj.top.filename, str(n_frames)))
    state.run()
    return state


t0 = traj[:]

iter_(t0)
state = rms2d_()

# make sure equal result: YES
