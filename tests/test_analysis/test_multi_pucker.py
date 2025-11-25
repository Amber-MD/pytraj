import pytraj as pt
from pytraj.testing import aa_eq, tempfolder
from utils import tz2_ortho_trajin, tz2_ortho_top


def test_multi_pucker():
    cm = """
    parm {}
    trajin {}
    multi_pucker :1-10 resrange 1-5 method altona range360 amplitude
    """.format(tz2_ortho_top, tz2_ortho_trajin)

    with tempfolder():
        # Run cpptraj state
        state = pt.datafiles.load_cpptraj_state(cm).run()
        cpptraj_results = state.data[:]

        # Run pytraj's multi_pucker
        traj = pt.iterload(tz2_ortho_trajin, tz2_ortho_top)
        pytraj_results = pt.multi_pucker(
            traj,
            mask=":1-10",
            resrange="1-5",
            method="altona",
            range360=True,
            amplitude=True,
        )

        # Compare results
        aa_eq(pytraj_results, cpptraj_results)