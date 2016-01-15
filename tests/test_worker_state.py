#!/usr/bin/env python
from __future__ import print_function
import unittest
from functools import partial
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.parallel.base import worker_state
from pytraj.tools import concat_dict


class TestWorkerState(unittest.TestCase):
    '''temp test. name will be changed.
    '''

    def test_different_cores(self):
        # use different frame_slice with different n_cores to test
        REF_0 = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")[0]

        for frame_slice in [(0, 100),
                            (10, 100, 3),
                            (50, 80, 2),
                            (51, 80, 3), ]:
            traj = pt.iterload("./data/tz2.nc",
                               "./data/tz2.parm7",
                               frame_slice=frame_slice)
            saved_angle = pt.angle(traj, ':3 :7 :9')
            saved_dist = pt.distance(traj, ':2 :10')
            saved_rmsd = pt.rmsd(traj, ref=REF_0, mask='@CA')

            lines = ['angle :3 :7 :9', 'distance :2 :10']

            for n_cores in [2, 3]:
                data_list = [worker_state(rank, n_cores, traj, lines)
                             for rank in range(n_cores)]
                final_data = concat_dict([x[1] for x in data_list])
                aa_eq(final_data['Ang_00002'], saved_angle)
                aa_eq(final_data['Dis_00003'], saved_dist)

    def test_multiple_cores(self):
        from multiprocessing import Pool
        traj = pt.iterload('data/tz2.nc', 'data/tz2.parm7')
        for _ in range(10):
            traj._load(traj.filelist)
        saved_angle = pt.angle(traj, ':3 :10 :11')
        saved_dist = pt.distance(traj, ':3 :10')

        for n_cores in [2, 3]:
            lines = ['angle :3 :10 :11', 'distance :3 :10']
            pfuncs = partial(worker_state,
                             n_cores=n_cores,
                             traj=traj,
                             dtype='dict',
                             lines=lines)
            p = Pool(n_cores)
            data_list = p.map(pfuncs, [rank for rank in range(n_cores)])
            p.close()
            p.join()
            data_list_sorted_rank = (data[1]
                                     for data in sorted(data_list,
                                                        key=lambda x: x[0]))
            final_data = concat_dict(data_list_sorted_rank)
            aa_eq(final_data['Ang_00002'], saved_angle)
            aa_eq(final_data['Dis_00003'], saved_dist)


if __name__ == "__main__":
    unittest.main()
