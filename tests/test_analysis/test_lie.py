import sys
import pytraj as pt
from utils import fn
from pytraj.testing import cpptraj_test_dir, aa_eq


def test_lie():
    topology_fn = cpptraj_test_dir + '/FtuFabI.NAD.TCL.parm7'
    trajin_fn = cpptraj_test_dir + '/FtuFabI.NAD.TCL.nc'
    traj = pt.iterload(trajin_fn, topology_fn)
    cpptraj_in = """
    parm {}
    trajin {}
    lie LIE :TCS cutvdw 12 cutelec 12
    """.format(topology_fn, trajin_fn)
    state = pt.load_cpptraj_state(cpptraj_in)
    state.run()

    mask = ':TCS'
    options = 'cutvdw 12 cutelec 12'
    data = pt.lie(traj, mask, options=options)

    for key in ['LIE[EELEC]', 'LIE[EVDW]']:
        aa_eq(data[key], state.data[key])

    if sys.platform.startswith('linux'):
        data_parallel = pt.pmap(pt.lie, traj, mask, options, n_cores=3)
        for key in ['LIE[EELEC]', 'LIE[EVDW]']:
            aa_eq(data_parallel[key], state.data[key])
