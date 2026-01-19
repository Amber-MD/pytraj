import unittest
from pytraj import *
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
import numpy as np


# tests to see if function can find an atom if we search with a coordinate
# very close to the atom's xyz coordinates
def test_search_for_atom():
    tz2_traj = pt.datafiles.load_tz2()

    # find coordinates of a given atom, make sure count_in_voxel finds
    # that atom at the right frame
    idx = 0
    frame = 0
    xyz = tz2_traj[frame].atom(idx)
    pop = pt.count_in_voxel(tz2_traj, "", xyz, 3)
    assert (idx in pop[frame])

    idx = 3
    frame = 4
    xyz = tz2_traj[frame].atom(idx)
    pop = pt.count_in_voxel(tz2_traj, "", xyz, 3)
    assert (idx in pop[frame])

# tests to make sure function doesn't put an atom in the list if its
# not in the voxel
def test_voxel_doesnt_contain():
    tz2_traj = pt.datafiles.load_tz2()
    idx = 0
    x, y, z = tz2_traj[0].atom(idx)
    size = 3
    voxel = (x + size, y + size, z + size)
    pop = pt.count_in_voxel(tz2_traj, "", voxel, size)

    assert (idx not in pop[0])

# example gist analysis for voxel centered at (35.26, 38.23, 1.66) with edge length 10
def test_gist_survival_ex():
    tz2_traj = pt.datafiles.load_tz2_ortho()
    wat_atoms = tz2_traj.top.select(":WAT")
    pop = pt.count_in_voxel(tz2_traj, ":WAT@O",
                            (35.26, 38.23, 1.66), 10)
    # print water molecules in voxel at frame 0.
    #
    # For survival time correlation function,
    # plot the number of these particular beginning frame water atoms that remain in the voxel
    # throughout the course of the simulation (if the voxel is part of the bulk solvent,
    # the function will decay very quickly, but if it is near the interface, then decay
    # will be slower as waters will stay in enthalpically favorable points).
    orig_frame_waters = set(pop[0])
    # NOTE: Reader will need to modify this to also exclude waters that leave in an intermediate frame
    # but then return later. (e.g. water 250 is in frame 1, but not in frame 2, then back in frame 3
    # it does not count as a "surviving" water that has remained in the voxel.
    survival_decay = [
        len(orig_frame_waters.intersection(set(pop[i])))
        for i in range(len(pop))
    ]


def test_count_in_voxel_frame_indices():
    """Test that frame_indices parameter works correctly for count_in_voxel"""
    tz2_traj = pt.datafiles.load_tz2()

    # Test with specific frame indices
    frame_indices = [0, 2, 4]
    voxel_center = (0, 0, 0)
    voxel_size = 10

    # Get results with frame_indices
    result_with_indices = pt.count_in_voxel(tz2_traj, '@CA', voxel_center, voxel_size, frame_indices=frame_indices)

    # Get results without frame_indices (all frames)
    result_all_frames = pt.count_in_voxel(tz2_traj, '@CA', voxel_center, voxel_size)

    # Check that we got the right number of results
    assert len(result_with_indices) == len(frame_indices), f"Expected {len(frame_indices)} results, got {len(result_with_indices)}"
    assert len(result_all_frames) == len(tz2_traj), f"Expected {len(tz2_traj)} results, got {len(result_all_frames)}"

    # Check that results for specific frames match
    for i, frame_idx in enumerate(frame_indices):
        assert result_with_indices[i] == result_all_frames[frame_idx], f"Frame {frame_idx} results don't match"


