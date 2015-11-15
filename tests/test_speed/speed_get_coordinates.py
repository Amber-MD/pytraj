#!/usr/bin/env python
import pytraj as pt
from pytraj.testing import Timer, aa_eq

n_frames = 5000
traj = pt.iterload('md.nc', 'tc5bwat.top', frame_slice=(0, n_frames))


@Timer()
def _get_coordinates(traj):
    return pt.get_coordinates(traj,
                              autoimage=True,
                              rmsfit=(0, '@CA'),
                              mask='@CA')


@Timer()
def run_cpptraj():
    text = '''
    parm tc5bwat.top
    trajin md.nc 0 %s
    autoimage
    rms first @CA
    strip !@CA
    createcrd
    ''' % str(n_frames)
    state = pt.datafiles.load_cpptraj_state(text)
    state.run()
    state.data[-1].xyz


_get_coordinates(traj)
run_cpptraj()

# assert True? Yes
# result: pytraj 11 s, cpptraj 8 s for 5000 frames
