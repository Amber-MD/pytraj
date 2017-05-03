import pytraj as pt
from pytraj.testing import aa_eq

# local
from utils import fn


def test_multidihedral():
    trajin = fn('Tc5b.x')
    topin = fn('Tc5b.top')
    traj = pt.iterload(trajin, top=topin)
    command = "resrange 2-19 phi psi"
    state = pt.load_cpptraj_state('''
    parm {}
    trajin {}
    multidihedral resrange 2-19 phi psi
    '''.format(topin, trajin))
    state.run()
    cpp_out = state.data[1:].to_dict()
    out = pt.multidihedral(traj, command, dtype='dict')
    for key, value in out.items():
        aa_eq(value, cpp_out.get(key))
