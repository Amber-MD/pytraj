import pytraj as pt
from pytraj.testing import tempfolder, assert_equal_dict
from utils import tz2_ortho_trajin, tz2_ortho_top


def test_tordiff():
    cm = f"""
    parm {tz2_ortho_top}
    trajin {tz2_ortho_trajin}
    tordiff TOR :WAT@O out tor.dat
    """

    with tempfolder():
        state = pt.datafiles.load_cpptraj_state(cm).run()
        cpptraj_results = state.data.to_dict()

        traj = pt.iterload(tz2_ortho_trajin, tz2_ortho_top)
        pytraj_results = pt.tordiff(
            traj,
            mask=":WAT@O",
            out="tor.dat"
        )

        cpptraj_results.pop('tz2.ortho.parm7', None)
        assert_equal_dict(pytraj_results, cpptraj_results)