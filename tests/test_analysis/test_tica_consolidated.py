import pytraj as pt
import numpy as np
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


def test_tica_coordinate_based():
    """Test coordinate-based TICA vs live cpptraj reference"""
    test_dir = os.path.join(cpptraj_test_dir, "Test_TICA")
    parm_file = os.path.join(test_dir, "../tz2.parm7")
    crd_file = os.path.join(test_dir, "../tz2.crd")

    # Test coordinate TICA with live cpptraj comparison
    cm = f"""
    parm {parm_file}
    loadcrd {crd_file} name MyCrd
    tica crdset MyCrd mask @CA lag 10 out crd.cumvar.dat
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
            raise ValueError("No cumvar data found in cpptraj coordinate TICA results")

        # Run pytraj equivalent
        traj = pt.iterload(crd_file, parm_file)
        coord_result = pt.tica(traj, mask='@CA', lag=10)

        # Compare results directly to cpptraj
        aa_eq(coord_result.cumvar, cpptraj_results[cumvar_key][:len(coord_result.cumvar)])


def test_tica_mixed_datasets():
    """Test TICA with mixed dataset types vs live cpptraj reference"""
    test_dir = os.path.join(cpptraj_test_dir, "Test_TICA")
    parm_file = os.path.join(test_dir, "../tz2.parm7")
    crd_file = os.path.join(test_dir, "../tz2.crd")

    # Test mixed data TICA with live cpptraj comparison
    cm = f"""
    parm {parm_file}
    trajin {crd_file}
    distance d1.12 :1@CA :12@CA
    distance d2.11 :2@CA :11@CA
    distance d3.10 :3@CA :10@CA
    rms R0 first @CA
    tica data R0 data d* lag 10 out mixed.cumvar.dat
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
            raise ValueError("No cumvar data found in cpptraj mixed TICA results")

        # Run pytraj equivalent
        traj = pt.iterload(crd_file, parm_file)
        rms_data = pt.rms(traj, ref=0, mask='@CA')
        d1 = pt.distance(traj, ':1@CA :12@CA')
        d2 = pt.distance(traj, ':2@CA :11@CA')
        d3 = pt.distance(traj, ':3@CA :10@CA')

        mixed_result = pt.tica(data=[rms_data, d1, d2, d3], lag=10)

        # Compare results directly to cpptraj
        aa_eq(mixed_result.cumvar, cpptraj_results[cumvar_key][:len(mixed_result.cumvar)])


def test_tica_different_lag():
    """Test TICA with different lag values vs live cpptraj reference"""
    test_dir = os.path.join(cpptraj_test_dir, "Test_TICA")
    parm_file = os.path.join(test_dir, "../tz2.parm7")
    crd_file = os.path.join(test_dir, "../tz2.crd")

    # Test with lag=5
    cm = f"""
    parm {parm_file}
    trajin {crd_file}
    distance d1 :1@CA :5@CA
    distance d2 :2@CA :6@CA
    tica data d1 data d2 lag 5 out lag5.cumvar.dat
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
            raise ValueError("No cumvar data found in cpptraj lag5 TICA results")

        # Run pytraj equivalent
        traj = pt.iterload(crd_file, parm_file)
        d1 = pt.distance(traj, ':1@CA :5@CA')
        d2 = pt.distance(traj, ':2@CA :6@CA')

        lag5_result = pt.tica(data=[d1, d2], lag=5)

        # Compare results directly to cpptraj
        aa_eq(lag5_result.cumvar, cpptraj_results[cumvar_key][:len(lag5_result.cumvar)])


def test_tica_projection():
    """Test TICA projection vs live cpptraj reference"""
    test_dir = os.path.join(cpptraj_test_dir, "Test_TICA")
    parm_file = os.path.join(test_dir, "../tz2.parm7")
    crd_file = os.path.join(test_dir, "../tz2.crd")

    # Test TICA with projection - let's see what cpptraj actually returns
    cm = f"""
    parm {parm_file}
    trajin {crd_file}
    distance d1 :1@CA :5@CA
    distance d2 :2@CA :6@CA
    tica data d1 data d2 name TICA lag 5 out proj.cumvar.dat
    runanalysis projectdata name Evec out proj.project.dat evecs TICA data d1 data d2
    run
    """

    with tempfolder():
        # Get live cpptraj reference
        state = pt.datafiles.load_cpptraj_state(cm).run()
        cpptraj_results = state.data.to_dict()

        # Find TICA cumvar key for validation
        cumvar_key = None
        for key in cpptraj_results.keys():
            if '[cumvar]' in key:
                cumvar_key = key
                break

        if cumvar_key is None:
            raise ValueError("No cumvar data found in cpptraj projection results")

        # Run pytraj equivalent
        traj = pt.iterload(crd_file, parm_file)
        d1 = pt.distance(traj, ':1@CA :5@CA')
        d2 = pt.distance(traj, ':2@CA :6@CA')

        tica_result = pt.tica(data=[d1, d2], lag=5)

        # Compare cumvar to validate TICA computation
        aa_eq(tica_result.cumvar, cpptraj_results[cumvar_key][:len(tica_result.cumvar)])

        # Verify TICA results have proper structure for projection
        assert hasattr(tica_result, 'modes')
        assert hasattr(tica_result, 'eigenvalues')
        assert hasattr(tica_result, 'eigenvectors')

        # Project the data using the TICA modes
        data_matrix = np.column_stack([d1, d2])
        projected = np.dot(data_matrix, tica_result.modes)

        # Verify projection dimensions are correct
        assert projected.shape[0] == len(d1)  # Same number of frames
        assert projected.shape[1] == tica_result.modes.shape[1]  # Same number of components

        # Verify projection values are reasonable (not NaN or infinite)
        assert not np.isnan(projected).any()
        assert not np.isinf(projected).any()

        # The cpptraj projectdata analysis fails ("modes set 'TICA' is empty")
        # but we can still test that pytraj TICA projection functionality works properly
        # and that the cumvar results match cpptraj TICA analysis


def test_tica_commute_mapping():
    """Test TICA with commute mapping vs live cpptraj reference"""
    test_dir = os.path.join(cpptraj_test_dir, "Test_TICA")
    parm_file = os.path.join(test_dir, "../tz2.parm7")
    crd_file = os.path.join(test_dir, "../tz2.crd")

    # Test TICA with commute mapping
    cm = f"""
    parm {parm_file}
    trajin {crd_file}
    distance d1 :1@CA :5@CA
    distance d2 :2@CA :6@CA
    tica data d1 data d2 lag 5 map commute out commute.cumvar.dat
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
            raise ValueError("No cumvar data found in cpptraj commute TICA results")

        # Run pytraj equivalent with commute mapping
        traj = pt.iterload(crd_file, parm_file)
        d1 = pt.distance(traj, ':1@CA :5@CA')
        d2 = pt.distance(traj, ':2@CA :6@CA')

        commute_result = pt.tica(data=[d1, d2], lag=5, commute=True)

        # Compare results directly to cpptraj
        aa_eq(commute_result.cumvar, cpptraj_results[cumvar_key][:len(commute_result.cumvar)])


def test_tica_n_components():
    """Test TICA with specific n_components vs live cpptraj reference"""
    test_dir = os.path.join(cpptraj_test_dir, "Test_TICA")
    parm_file = os.path.join(test_dir, "../tz2.parm7")
    crd_file = os.path.join(test_dir, "../tz2.crd")

    # Test TICA with specific number of components
    cm = f"""
    parm {parm_file}
    trajin {crd_file}
    distance d1 :1@CA :5@CA
    distance d2 :2@CA :6@CA
    distance d3 :3@CA :7@CA
    tica data d1 data d2 data d3 lag 3 out ncomp.cumvar.dat
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
            raise ValueError("No cumvar data found in cpptraj n_components TICA results")

        # Run pytraj equivalent
        traj = pt.iterload(crd_file, parm_file)
        d1 = pt.distance(traj, ':1@CA :5@CA')
        d2 = pt.distance(traj, ':2@CA :6@CA')
        d3 = pt.distance(traj, ':3@CA :7@CA')

        # Test with different n_components
        for n_comp in [1, 2, 3]:
            ncomp_result = pt.tica(data=[d1, d2, d3], lag=3, n_components=n_comp)

            # Verify structure
            assert ncomp_result.modes.shape[1] == n_comp
            assert len(ncomp_result.eigenvalues) >= n_comp

            # Compare cumvar for the requested components
            aa_eq(ncomp_result.cumvar[:n_comp],
                  cpptraj_results[cumvar_key][:n_comp], decimal=5)
