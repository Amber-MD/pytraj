#!/usr/bin/env python
import pytraj as pt
from pytraj.testing import Timer, aa_eq

n_frames = 1000
traj = pt.iterload('GAAC3.5000frames.nc', 'GAAC.topo', frame_slice=(0, n_frames))

mask = '!@H='

@Timer()
def _get_coordinates(traj):
    return pt.get_coordinates(traj, autoimage=True, rmsfit=(0, mask), mask=mask)

@Timer()
def run_cpptraj():
    text = '''
    parm GAAC.topo
    trajin GAAC3.5000frames.nc 0 {0}
    autoimage
    rms first {1}
    strip !@H=
    createcrd
    '''.format(str(n_frames), mask)
    state = pt.datafiles.load_cpptraj_state(text)
    state.run()
    state.data[-1].xyz

_get_coordinates(traj)
run_cpptraj()

# assert True? Yes
#saved_xyz = pt.iterload('test.nc', 'test.tc5bwat.top').xyz
#aa_eq(saved_xyz, xyz)

# result: pytraj 11 s, cpptraj 8 s for 5000 frames
