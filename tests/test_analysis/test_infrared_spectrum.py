import pytraj as pt
from pytraj.testing import aa_eq, tempfolder, cpptraj_test_dir

def test_infraredspec():
    parm_file = f"{cpptraj_test_dir}/Test_systemVF/systemVF.parm7"
    traj_file = f"{cpptraj_test_dir}/Test_systemVF/systemVF.nc"

    cm = f"""
    parm {parm_file}
    trajin {traj_file}
    infraredspec IR out irspec.dat maxlag 5 tstep 0.1 rawout raw.dat
    """
    with tempfolder():
        # Run cpptraj state
        state = pt.datafiles.load_cpptraj_state(cm).run()
        cpptraj_results = state.data.to_dict()

        # Run pytraj's infraredspec
        traj = pt.iterload(traj_file, parm_file)
        pytraj_results = pt.infraredspec(
            traj,
            mask="IR",
            out="irspec.dat",
            maxlag=5,
            tstep=0.1,
            rawout="raw.dat",
        )

        # Compare results
        for k, v in pytraj_results.items():
            aa_eq(v, cpptraj_results[k])