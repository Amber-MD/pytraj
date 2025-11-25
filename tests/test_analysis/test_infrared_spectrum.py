import pytraj as pt
from pytraj.testing import aa_eq, tempfolder, cpptraj_test_dir

def test_infrared_spectrum():
    parm_file = f"{cpptraj_test_dir}/Test_systemVF/systemVF.parm7"
    traj_file = f"{cpptraj_test_dir}/Test_systemVF/systemVF.nc"

    cm = f"""
    parm {parm_file}
    trajin {traj_file}
    infraredspec IR out irspec.dat maxlag 5 tstep 0.1 rawout raw.dat
    """
    pt._verbose()
    with tempfolder():
        # Run cpptraj state
        state = pt.datafiles.load_cpptraj_state(cm).run()
        cpptraj_results = state.data[:]

        # Run pytraj's infrared_spectrum
        traj = pt.iterload(traj_file, parm_file)
        pytraj_results = pt.infrared_spectrum(
            traj,
            mask="IR",
            out="irspec.dat",
            maxlag=5,
            tstep=0.1,
            rawout="raw.dat",
        )

        # Compare results
        aa_eq(pytraj_results, cpptraj_results)