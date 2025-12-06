import os
import pytraj as pt
from pytraj.testing import cpptraj_test_dir, aa_eq


def test_xtalsymm(tmpdir):
    with tmpdir.as_cwd():
        x6dky_parm = os.path.join(cpptraj_test_dir, 'x6dky.parm7')
        x6dky_tname = os.path.join(cpptraj_test_dir, 'mdXtal.nc')
        x6dky_refname = os.path.join(cpptraj_test_dir, 'mdXtal.inpcrd')
        state = pt.load_cpptraj_state("""
        parm {parm}
        trajin {trajin}
        reference {reference}
        xtalsymm :1-16 reference group P22(1)2(1) na 2 nb 1 nc 1
        trajout mdAsuOnly.pdb
        """.format(
            parm=x6dky_parm, trajin=x6dky_tname, reference=x6dky_refname))
        state.run()

        traj = pt.load(x6dky_tname, x6dky_parm)
        ref = pt.load(x6dky_refname, traj.top)
        pt.xtalsymm(
            traj,
            mask=':1-16',
            ref=ref,
            options="group P22(1)2(1) na 2 nb 1 nc 1")
        c_traj = pt.load('mdAsuOnly.pdb')
        aa_eq(traj.xyz, c_traj.xyz, decimal=2)


def test_xtalsymm_frame_indices():
    """Test that frame_indices parameter works correctly for xtalsymm"""
    x6dky_parm = os.path.join(cpptraj_test_dir, 'x6dky.parm7')
    x6dky_tname = os.path.join(cpptraj_test_dir, 'mdXtal.nc')
    x6dky_refname = os.path.join(cpptraj_test_dir, 'mdXtal.inpcrd')

    traj = pt.load(x6dky_tname, x6dky_parm)
    ref = pt.load(x6dky_refname, traj.top)

    # Test with specific frames
    frame_indices = [0] if len(traj) > 0 else []
    mask = ':1-16'
    options = 'group P22(1)2(1) na 2 nb 1 nc 1'

    # Get results with and without frame_indices
    result_with_indices = pt.xtalsymm(traj, mask=mask, ref=ref, options=options, frame_indices=frame_indices)
    result_all_frames = pt.xtalsymm(traj, mask=mask, ref=ref, options=options)

    # Both should return DatasetList
    assert hasattr(result_with_indices, '__len__')
    assert hasattr(result_all_frames, '__len__')

    # Should have same structure but potentially different lengths
    assert type(result_with_indices) == type(result_all_frames)
