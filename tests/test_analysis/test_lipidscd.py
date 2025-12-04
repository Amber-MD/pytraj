import pytraj as pt
from pytraj.testing import aa_eq
import numpy as np

# local
from utils import fn


def test_lipidscd():
    parm = fn('DOPC.parm7')
    trajin = fn('DOPC.rst7')

    mask = ':OL,PC'
    cm = """
    parm {}
    trajin {}
    lipidscd {}
    """.format(parm, trajin, mask)

    state = pt.load_cpptraj_state(cm)
    state.run()
    expected_data = state.data[1:].to_dict()

    traj = pt.load(trajin, parm)
    out = pt.lipidscd(traj, mask)
    for k, v in out.items():
        aa_eq(v, expected_data[k])


def test_lipidscd_frame_indices():
    """Test that frame_indices parameter works correctly for lipidscd"""
    parm = fn('DOPC.parm7')
    trajin = fn('DOPC.rst7')

    traj = pt.load(trajin, parm)
    mask = ':OL,PC'

    # Test with specific frames
    frame_indices = [0] if len(traj) > 0 else []

    # Get results with and without frame_indices
    result_with_indices = pt.lipidscd(traj, mask, frame_indices=frame_indices)
    result_all_frames = pt.lipidscd(traj, mask)

    # Verify we get dict results
    assert isinstance(result_with_indices, dict)
    assert isinstance(result_all_frames, dict)

    # Both should have same keys
    assert result_with_indices.keys() == result_all_frames.keys()

    # Results with frame_indices should have expected length
    for key in result_with_indices.keys():
        assert len(result_with_indices[key]) == len(frame_indices)
        assert len(result_all_frames[key]) == len(traj)
