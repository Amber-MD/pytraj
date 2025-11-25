import pytraj as pt
from pytraj.testing import aa_eq, tempfolder
from utils import tc5b_trajin, tc5b_top

def test_toroidal_diffusion():
    cm = """
    parm {}
    trajin {}
    toroidal_diffusion :1-10 out test_diff_avg.dat diffout test_diff_results.dat time 2.0
    """.format(tc5b_top, tc5b_trajin)

    with tempfolder():
        # Run cpptraj state
        state = pt.datafiles.load_cpptraj_state(cm).run()
        cpptraj_results = state.data[:]

        # Run pytraj's toroidal_diffusion
        traj = pt.iterload(tc5b_trajin, tc5b_top)
        pytraj_results = pt.toroidal_diffusion(
            traj,
            mask=":1-10",
            out="test_diff_avg.dat",
            diffout="test_diff_results.dat",
            time=2.0,
        )

        # Compare results
        aa_eq(pytraj_results, cpptraj_results)