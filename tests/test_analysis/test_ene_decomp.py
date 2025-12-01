import pytraj as pt
from pytraj.testing import aa_eq, tempfolder, cpptraj_test_dir


def test_ene_decomp():
    # Use cpptraj_test_dir to locate the test files
    parm_file = f"{cpptraj_test_dir}/Test_Ewald/nacl.box.parm7"
    traj_file = f"{cpptraj_test_dir}/Test_Ewald/nacl.box.rst7"

    # Test 1: Basic energy decomposition
    cm = f"""
    parm {parm_file}
    trajin {traj_file}
    enedecomp ATM * out decomp.nacl.box.dat
    """

    with tempfolder():
        # Run cpptraj state
        state = pt.datafiles.load_cpptraj_state(cm).run()
        cpptraj_results = state.data[-1].values

        # Run pytraj's ene_decomp
        traj = pt.iterload(traj_file, parm_file)
        pytraj_results = pt.ene_decomp(
            traj,
            mask="*",
        )[-1]

        # Compare results
        aa_eq(pytraj_results, cpptraj_results)

    # Test 2: Energy decomposition with savecomponents
    cm_savecomponents = f"""
    parm {parm_file}
    trajin {traj_file}
    enedecomp ATN * out decomp.components.nacl.box.dat savecomponents
    """

    with tempfolder():
        # Run cpptraj state
        state = pt.datafiles.load_cpptraj_state(cm_savecomponents).run()
        cpptraj_results_savecomponents = state.data[-1].values

        # Run pytraj's ene_decomp with savecomponents
        traj = pt.iterload(traj_file, parm_file)
        pytraj_results_savecomponents = pt.ene_decomp(
            traj,
            mask="*",
            savecomponents=True,
        )[-1]

        # Compare results
        aa_eq(pytraj_results_savecomponents, cpptraj_results_savecomponents)

    # Test 3: Energy decomposition with a specific atom mask
    cm_mask = f"""
    parm {parm_file}
    trajin {traj_file}
    enedecomp ATM @4 out decomp.mask.nacl.box.dat
    """

    with tempfolder():
        # Run cpptraj state
        state = pt.datafiles.load_cpptraj_state(cm_mask).run()
        cpptraj_results_mask = state.data[-1].values

        # Run pytraj's ene_decomp with a specific mask
        traj = pt.iterload(traj_file, parm_file)
        pytraj_results_mask = pt.ene_decomp(
            traj,
            mask="@4",
        )[-1]

        # Compare results
        aa_eq(pytraj_results_mask, cpptraj_results_mask)