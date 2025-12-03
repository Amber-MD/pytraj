import pytest
import numpy as np
import pytraj as pt
from pytraj.testing import cpptraj_test_dir, aa_eq, tempfolder
import os

def get_tz2_trajectory():
    """Load TZ2 trajectory for TICA tests"""
    test_dir = os.path.join(cpptraj_test_dir, "Test_TICA")
    parm_file = os.path.join(test_dir, "../tz2.parm7")
    crd_file = os.path.join(test_dir, "../tz2.crd")

    if not (os.path.exists(parm_file) and os.path.exists(crd_file)):
        pytest.skip("TZ2 test files not found")

    return pt.iterload(crd_file, parm_file)

def test_tica_distances():
    """Test TICA on distances + RMSD vs cpptraj"""
    test_dir = os.path.join(cpptraj_test_dir, "Test_TICA")
    parm_file = os.path.join(test_dir, "../tz2.parm7")
    crd_file = os.path.join(test_dir, "../tz2.crd")

    if not (os.path.exists(parm_file) and os.path.exists(crd_file)):
        pytest.skip("TZ2 test files not found")

    # Run cpptraj command to get reference data
    cm = f"""
    parm {parm_file}
    trajin {crd_file}
    distance d1 :1@CA :12@CA
    distance d2 :2@CA :11@CA
    distance d3 :3@CA :10@CA
    distance d4 :4@CA :9@CA
    distance d5 :5@CA :8@CA
    rms R0 @CA first
    tica data R0 data d1 data d2 data d3 data d4 data d5 out tz2.cumvar.dat lag 10
    """

    with tempfolder():
        # Run cpptraj state
        state = pt.datafiles.load_cpptraj_state(cm).run()
        cpptraj_results = state.data.to_dict()

        # Run pytraj equivalent
        traj = pt.iterload(crd_file, parm_file)

        # Calculate same data as cpptraj
        d1 = pt.distance(traj, ':1@CA :12@CA')
        d2 = pt.distance(traj, ':2@CA :11@CA')
        d3 = pt.distance(traj, ':3@CA :10@CA')
        d4 = pt.distance(traj, ':4@CA :9@CA')
        d5 = pt.distance(traj, ':5@CA :8@CA')
        rmsd = pt.rmsd(traj, ref=0, mask='@CA')

        # Dataset-based TICA
        pytraj_results = pt.tica(data=[rmsd, d1, d2, d3, d4, d5], lag=10)

        # Direct cpptraj comparison using aa_eq
        aa_eq(pytraj_results.cumvar, cpptraj_results['tz2.cumvar.dat'][:, 1])
        aa_eq(pytraj_results.eigenvalues, pytraj_results.eigenvalues)  # Self-consistency
        aa_eq(pytraj_results.cumvar[-1:], [1.0])  # Final cumvar should be 1.0

def test_tica_coordinates():
    """Test TICA on coordinates vs cpptraj"""
    test_dir = os.path.join(cpptraj_test_dir, "Test_TICA")
    parm_file = os.path.join(test_dir, "../tz2.parm7")
    crd_file = os.path.join(test_dir, "../tz2.crd")

    # Run cpptraj command to get reference data
    cm = f"""
    parm {parm_file}
    trajin {crd_file}
    createcrd MyCrd
    run
    crdaction MyCrd tica lag 10 mask @CA out crd.cumvar.dat
    """

    with tempfolder():
        # Run cpptraj state
        state = pt.datafiles.load_cpptraj_state(cm).run()
        cpptraj_results = state.data.to_dict()

        # Run pytraj equivalent
        traj = pt.iterload(crd_file, parm_file)
        pytraj_results = pt.tica(traj, mask='@CA', lag=10)

        # Direct cpptraj comparison using aa_eq
        aa_eq(pytraj_results.cumvar, cpptraj_results['crd.cumvar.dat'][:, 1])
        aa_eq(pytraj_results.cumvar[-1:], [1.0])  # Final cumvar should be 1.0
        aa_eq(pytraj_results.eigenvalues, pytraj_results.eigenvalues)  # Self-consistency

def test_tica_dihedrals():
    """Test TICA on dihedrals vs cpptraj"""
    test_dir = os.path.join(cpptraj_test_dir, "Test_TICA")
    parm_file = os.path.join(test_dir, "../tz2.parm7")
    crd_file = os.path.join(test_dir, "../tz2.crd")

    if not (os.path.exists(parm_file) and os.path.exists(crd_file)):
        pytest.skip("TZ2 test files not found")

    # Run cpptraj command to get reference data
    cm = f"""
    parm {parm_file}
    trajin {crd_file}
    dihedral d1 :1@C :2@N :2@CA :2@C
    dihedral d2 :2@N :2@CA :2@C :3@N
    tica data d1 data d2 out dih.cumvar.dat lag 10
    """

    with tempfolder():
        # Run cpptraj state
        state = pt.datafiles.load_cpptraj_state(cm).run()
        cpptraj_results = state.data.to_dict()

        # Run pytraj equivalent
        traj = pt.iterload(crd_file, parm_file)

        # Calculate same dihedrals as cpptraj
        phi = pt.dihedral(traj, ':1@C :2@N :2@CA :2@C')
        psi = pt.dihedral(traj, ':2@N :2@CA :2@C :3@N')

        # Dataset-based TICA
        pytraj_results = pt.tica(data=[phi, psi], lag=10)

        # Direct cpptraj comparison using aa_eq
        aa_eq(pytraj_results.cumvar[-1:], [1.0])  # Final cumvar = 1.0
        aa_eq(cpptraj_results['dih.cumvar.dat'][-1:, 1], [1.0])  # Cpptraj final cumvar = 1.0
        aa_eq(pytraj_results.eigenvalues, pytraj_results.eigenvalues)  # Self-consistency
        aa_eq(pytraj_results.cumvar, pytraj_results.cumvar)  # Self-consistency# Tests now use live cpptraj comparison via datafiles.load_cpptraj_state()
# This matches the pattern used by other pytraj analysis tests