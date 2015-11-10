import pytraj as pt
from pytraj.testing import Timer

traj = pt.iterload('GAAC3.1000frames.nc', 'GAAC.parm7', frame_slice=(0, 1000))
hb = pt.hbond(traj, '!:WAT')

masklist = hb._amber_mask()
cpp_cm = ['distance ' + _ for _ in masklist]
cpp_mask = '\n'.join(cpp_cm)

@Timer()
def test_pytraj():
    print('pytraj')
    traj = pt.iterload('GAAC3.1000frames.nc', 'GAAC.parm7', frame_slice=(0, 1000))
    data = pt.distance(traj, masklist)
    #print(data)

@Timer()
def test_state():
    print('cpptraj')
    state = pt.load_cpptraj_state('''
    parm GAAC.parm7
    trajin GAAC3.1000frames.nc 1 1000
    {0}'''.format(cpp_mask))
    state.run()
    #print(state.data[1:].values)

test_pytraj()
test_state()

# timing output (seconds) with 1000 frames.

# pytraj
# 4.944881916046143
# cpptraj
# 4.855885982513428
