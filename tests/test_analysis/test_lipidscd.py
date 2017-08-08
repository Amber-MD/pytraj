import pytraj as pt
from pytraj.testing import aa_eq

# local
from utils import fn


def test_lipidscd():
    parm = fn('DOPC.parm7')
    trajin = fn('DOPC.rst7')

    mask = ':OL,PC'
    cm = """
    parm {}
    trajin {}
    lipidscd {}
    """.format(parm, trajin, mask)

    state = pt.load_cpptraj_state(cm)
    state.run()
    expected_data = state.data[1:].to_dict()

    traj = pt.load(trajin, parm)
    out = pt.lipidscd(traj, mask)
    for k, v in out.items():
        aa_eq(v, expected_data[k])
