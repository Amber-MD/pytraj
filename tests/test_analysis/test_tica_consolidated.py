import pytest
import pytraj as pt
from pytraj.testing import cpptraj_test_dir, aa_eq
import os
from pytraj.testing import tempfolder

def test_tica_simple():
    """Test TICA vs live cpptraj reference data"""
    test_dir = os.path.join(cpptraj_test_dir, "Test_TICA")
    parm_file = os.path.join(test_dir, "../tz2.parm7")
    crd_file = os.path.join(test_dir, "../tz2.crd")

    # Test dihedral TICA with live cpptraj comparison
    cm = f"""
    parm {parm_file}
    trajin {crd_file}
    dihedral d1 :1@C :2@N :2@CA :2@C
    dihedral d2 :2@N :2@CA :2@C :3@N
    tica data d1 data d2 out dih.cumvar.dat lag 10
    """

    with tempfolder():
        # Get live cpptraj reference
        state = pt.datafiles.load_cpptraj_state(cm).run()
        cpptraj_results = state.data.to_dict()

        # Find TICA cumvar key
        cumvar_key = None
        for key in cpptraj_results.keys():
            if '[cumvar]' in key:
                cumvar_key = key
                break

        if cumvar_key is None:
            raise ValueError("No cumvar data found in cpptraj results")

        # Run pytraj equivalent
        traj = pt.iterload(crd_file, parm_file)
        phi = pt.dihedral(traj, ':1@C :2@N :2@CA :2@C')
        psi = pt.dihedral(traj, ':2@N :2@CA :2@C :3@N')
        dih_result = pt.tica(data=[phi, psi], lag=10)

        # Compare results directly to cpptraj
        aa_eq(dih_result.cumvar, cpptraj_results[cumvar_key][:len(dih_result.cumvar)])
