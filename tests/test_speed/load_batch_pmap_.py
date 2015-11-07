#!/usr/bin/env python

import pytraj as pt
import os
from pytraj.testing import Timer

NCORES = int(os.environ.get('NCORES', 4))
print('NCORES = ', NCORES)

traj = pt.iterload('md.nc', 'tc5bwat.top')
traj.load(traj.filename)
traj.load(traj.filename)
traj.load(traj.filename)
traj.load(traj.filename)


@Timer()
def test_cpptraj(traj=traj):
    state = pt.load_batch(traj, '''
    autoimage
    distance :1 :3
    distance :5 :18
    molsurf @CA
    multidihedral
    ''')
    state.run()
    return state


@Timer()
def test_pmap(traj=traj, n_cores=NCORES):
    return pt._load_batch_pmap(
        n_cores=n_cores,
        traj=traj,
        lines=['autoimage', 'distance :1 :3', 'distance :5 :18', 'molsurf @CA',
               'multidihedral'])


test_pmap()
test_cpptraj()
