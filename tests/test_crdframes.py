#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestCrdFrames(unittest.TestCase):

    def test_crdframes(self):
        '''test crdframes in cpptraj
        '''
        max_frames = 50
        traj = pt.iterload('data/tz2.nc',
                           'data/tz2.parm7',
                           frame_slice=(0, max_frames, 2))

        state = pt.load_cpptraj_state('''
                parm {0}
                trajin {1} 1 {2} 2
                trajin {1} 1 {2} 2
                loadtraj name traj
                crdaction traj rms crdframes 1,10
                crdaction traj rms crdframes 1,30,2
                '''.format(traj.top.filename, traj.filename, max_frames))
        state.run()

        rmsd_0 = pt.rmsd(traj, ref=0, frame_indices=range(10))
        rmsd_crdframes = state.data[2].values
        aa_eq(rmsd_0, rmsd_crdframes)

        traj2 = traj.copy()
        assert traj2.n_frames == traj.n_frames, 'must have the same n_frames'
        traj2._load(traj.filename, frame_slice=(0, max_frames, 2))
        assert traj2.n_frames == 2 * traj.n_frames, 'n_frames must be doubled after reload'

        rmsd_1 = pt.rmsd(traj2, ref=0, frame_indices=range(0, 30, 2))
        rmsd_crdframes = state.data[3].values
        aa_eq(rmsd_1, rmsd_crdframes)


if __name__ == "__main__":
    unittest.main()
