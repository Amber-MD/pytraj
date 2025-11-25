import pytraj as pt
from pytraj.testing import aa_eq, tempfolder
from utils import tz2_ortho_trajin, tz2_ortho_top


def test_time_analysis():
    cm = """
    parm {}
    trajin {}
    time_analysis :1-10
    """.format(tz2_ortho_top, tz2_ortho_trajin)

    with tempfolder():
        # Run cpptraj state
        state = pt.datafiles.load_cpptraj_state(cm).run()
        cpptraj_results = state.data[:]

        # Run pytraj's time_analysis
        traj = pt.iterload(tz2_ortho_trajin, tz2_ortho_top)
        pytraj_results = pt.time_analysis(
            traj,
            mask=":1-10",
        )

        # Compare results
        aa_eq(pytraj_results, cpptraj_results)