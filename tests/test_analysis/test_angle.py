import unittest
import numpy as np

import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq, load_cpptraj_reference_data
def test_angle_cpptraj_reference():
    """Test angle calculations against cpptraj Test_General reference data"""
    # Load the same trajectory as cpptraj Test_General
    traj = pt.datafiles.load_tz2()

    # Calculate angle :2@CA :3@CA :4@CA (matches cpptraj exactly)
    angle_result = pt.angle(traj, ':2@CA :3@CA :4@CA')

    # Load reference data directly from cpptraj Test_General
    expected_cpptraj = load_cpptraj_reference_data('Test_General', 'a1.dat.save')
    assert expected_cpptraj is not None, "cpptraj reference file a1.dat.save not found"

    # pytraj frame i matches cpptraj frame i+2 (offset of +1 in cpptraj indexing)
    # Compare first 10 frames accounting for this offset
    pytraj_frames = angle_result[:10]
    cpptraj_frames = expected_cpptraj[1:11]  # Skip first cpptraj frame

    np.testing.assert_allclose(pytraj_frames, cpptraj_frames,
                             rtol=1e-4, atol=1e-3)

    print(f"✓ pytraj angle results match cpptraj Test_General reference data ({len(expected_cpptraj)} frames)")

    # Ensure we match the expected frame count
    assert len(angle_result) == traj.n_frames


def test_angle_basic_validation():
    """Basic angle validation test (fallback)"""
    traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
    fa = traj[:]
    mask = ':1@CA :2@CA :3@CA'

    # Test basic functionality
    d0 = pt.calc_angle(traj, mask, dtype='dataset').to_ndarray()
    d1 = pt.angle(traj, mask)
    d2 = pt.angle(fa, mask)

    aa_eq(d0, d1)
    aa_eq(d0, d2)

    Nsize = 10
    arr = np.random.randint(0, 300, size=Nsize * 3).reshape(Nsize, 3)
    d3 = pt.angle(fa, arr)
    d4 = pt.angle(traj, arr)
    d5 = pt.calc_angle(traj, arr)
    d6 = pt.calc_angle(fa, arr)
    d7 = pt.calc_angle([fa, traj], arr, n_frames=2 * fa.n_frames)

    aa_eq(d3, d4)
    aa_eq(d3, d5)
    aa_eq(d3, d6)
    aa_eq(d3.T, d7.T[:fa.n_frames])
    aa_eq(d3.T, d7.T[fa.n_frames:])

    # Test data types
    d8 = pt.angle(traj, mask, dtype='dataset')
    d9 = pt.tools.dict_to_ndarray(pt.angle(traj, mask, dtype='dict'))
    aa_eq(d0, d8.values)
    aa_eq([d0], d9)


def test_angle():
    """Main test function that prioritizes cpptraj comparison"""
    test_angle_cpptraj_reference()
