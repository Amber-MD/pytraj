import pytraj as pt
from pytraj.testing import aa_eq, tempfolder
from utils import tz2_ortho_trajin, tz2_ortho_top


def test_create_reservoir():
    cm = """
    parm {}
    trajin {}
    create_reservoir :1-10
    """.format(tz2_ortho_top, tz2_ortho_trajin)

    with tempfolder():
        # Run cpptraj state
        state = pt.datafiles.load_cpptraj_state(cm).run()
        cpptraj_results = state.data[:]

        # Run pytraj's create_reservoir
        traj = pt.iterload(tz2_ortho_trajin, tz2_ortho_top)
        pytraj_results = pt.create_reservoir(
            traj,
            mask=":1-10",
        )

        # Compare results
        aa_eq(pytraj_results, cpptraj_results)