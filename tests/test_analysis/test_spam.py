import os
import pytraj as pt
from pytraj.testing import aa_eq, cpptraj_test_dir
from pytraj.utils import tempfolder

# local
from utils import tz2_truncoct_trajin
from utils import tz2_truncoct_top


def test_spam():
    # FIXME: assert?
    peaks_xyz = os.path.join(cpptraj_test_dir, 'Test_SPAM', 'peaks.xyz')
    command = """
    parm {}
    trajin {}
    autoimage
    spam {} name SPAM cut 12.0 info spam.info out spam.dat reorder \
         summary summary.dat
    trajout test.mdcrd onlyframes 1-2
    """.format(tz2_truncoct_top, tz2_truncoct_trajin, peaks_xyz)

    state = pt.load_cpptraj_state(command)
    with tempfolder():
        state.run()

    traj = pt.iterload(tz2_truncoct_trajin, tz2_truncoct_top)
    # cm = 'SPAM cut 12.0 info spam.info out spam.dat reorder summary summary.dat'
    cm = 'SPAM cut 12.0 reorder'
    with tempfolder():
        spam_out = pt.spam(traj, peak_file=peaks_xyz, command=cm)
