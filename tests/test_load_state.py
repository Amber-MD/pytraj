#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj import load_batch


class TestState(unittest.TestCase):
    def test_loading(self):
        for frame_slice in [(0, -1, 1), (0, 8, 2)]:
            traj = pt.iterload('data/tz2.nc', './data/tz2.parm7',
                               frame_slice=frame_slice)

            text = '''
            rms @CA
            radgyr @CA nomax
            '''

            s = load_batch(traj, text)
            s.run()

            dslist = s.datasetlist
            rmsd0 = pt.rmsd(traj, 0, '@CA')
            r0 = pt.radgyr(traj, '@CA')
            aa_eq(rmsd0, dslist[0])
            aa_eq(r0, dslist[1])
            assert len(dslist[0]) == traj.n_frames

    def test_raise_if_not_trajiter(self):
        traj = pt.iterload('data/tz2.nc', './data/tz2.parm7')
        t0 = traj[:]

        text = '''
        rms @CA
        radgyr @CA nomax
        '''

        self.assertRaises(ValueError, lambda: load_batch(t0, text))


if __name__ == "__main__":
    unittest.main()
