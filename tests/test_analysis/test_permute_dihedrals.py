import pytraj as pt
from utils import fn


def test_permute_dihedrals(tmpdir):
    with tmpdir.as_cwd():
        traj = pt.load(fn('tz2.rst7'), fn('tz2.parm7'))
        pt.all_actions.permute_dihedrals(
            traj, 'hey.nc', options="interval -120 phi psi")
        trajout = pt.load('hey.nc', traj.top)
        assert trajout.n_frames == 70
