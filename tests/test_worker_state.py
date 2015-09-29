#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.parallel import _worker_state
from pytraj.tools import concat_dict


class TestWorkerState(unittest.TestCase):
    '''temp test. name will be changed.
    '''
    def test_different_cores(self):
        # use different frame_slice with different n_cores to test
        for frame_slice in [(0, 100),
                            (10, 100, 3),
                            (50, 80, 2),]:
            traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7", frame_slice=frame_slice)
            saved_angle = pt.angle(traj, ':3 :7 :9')
            saved_dist= pt.distance(traj, ':2 :10')

            lines = ['angle :3 :7 :9', 'distance :2 :10']

            for n_cores in [2, 3, 4, 5, 6, 7, 8]:
                data_list = [_worker_state(n_cores, rank, traj, lines)
                        for rank in range(n_cores)
                        ]
                data_list_sorted_rank = (data[1] for data in sorted(data_list, key=lambda x : x[0]))
                final_data = concat_dict(data_list_sorted_rank)
                aa_eq(final_data['Ang_00001'], saved_angle)
                aa_eq(final_data['Dis_00002'], saved_dist)

if __name__ == "__main__":
    unittest.main()
