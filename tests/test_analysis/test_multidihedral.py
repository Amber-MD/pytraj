import pytraj as pt
from pytraj.testing import aa_eq

# local
from utils import fn


def test_multidihedral():
    trajin = fn('Tc5b.x')
    topin = fn('Tc5b.top')
    traj = pt.iterload(trajin, top=topin)
    command = "resrange 2-19 phi psi"
    state = pt.load_cpptraj_state('''
    parm {}
    trajin {}
    multidihedral resrange 2-19 phi psi
    '''.format(topin, trajin))
    state.run()
    cpp_out = state.data[1:].to_dict()
    out = pt.multidihedral(traj, command, dtype='dict')
    for key, value in out.items():
        aa_eq(value, cpp_out.get(key))


def test_multidihedral_frame_indices():
    """Test that frame_indices parameter works correctly for multidihedral"""
    trajin = fn('Tc5b.x')
    topin = fn('Tc5b.top')
    traj = pt.iterload(trajin, top=topin)

    # Test with specific frame indices
    frame_indices = [0, 2, 4]
    command = "phi psi"

    # Get results with frame_indices
    result_with_indices = pt.multidihedral(traj, command, frame_indices=frame_indices)

    # Get results without frame_indices (all frames)
    result_all_frames = pt.multidihedral(traj, command)

    # Check that we got the right number of frames
    for key in result_with_indices.keys():
        assert len(result_with_indices[key]) == len(frame_indices), f"Dataset {key}: expected {len(frame_indices)} values, got {len(result_with_indices[key])}"
        assert len(result_all_frames[key]) == len(traj), f"Dataset {key}: expected {len(traj)} values, got {len(result_all_frames[key])}"

        # Check that specific frame values match
        for i, frame_idx in enumerate(frame_indices):
            aa_eq(result_with_indices[key][i], result_all_frames[key][frame_idx])
