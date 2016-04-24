#!/usr/bin/env python

from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestNormalDistance(unittest.TestCase):

    def test_general(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        mask = ':1@CA :14@CB'
        d0 = pt.distance(traj, mask)
        d1 = pt.distance(traj, mask)
        d2 = pt.distance(fa, mask)

        aa_eq(d0, d1)
        aa_eq(d0, d2)

        Nsize = 12
        arr = np.random.randint(0, 300, size=Nsize * 2).reshape(Nsize, 2)
        d3 = pt.distance(fa, arr)
        d4 = pt.distance(traj, arr)
        d5 = pt.distance(traj, arr)
        d6 = pt.distance(fa, arr)
        d7 = pt.distance([fa, traj], arr, n_frames=2 * fa.n_frames)
        d8 = pt.distance([fa, traj], arr, n_frames=2 * fa.n_frames, dtype='dataset')
        aa_eq(d3, d4)
        aa_eq(d3, d5)
        aa_eq(d3, d6)
        aa_eq(d3.T, d7.T[:fa.n_frames])
        aa_eq(d3.T, d7.T[fa.n_frames:])
        aa_eq(d7, d8.values)

        # raise
        self.assertRaises(ValueError, lambda: pt.dihedrals(traj, [[0, 3, 2]]))

    def test_calculate_distance_without_specifying_n_frames(self):
        # TrajectoryIterator
        import numpy as np
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        arr = pt.distance(traj(stop=4), [0, 5])
        arr1 = pt.distance(traj(stop=4), [0, 5], n_frames=4)
        assert np.all(arr == arr1)

        arr2 = pt.distance(traj(stop=1000), [0, 5])
        arr3 = pt.distance(traj(stop=traj.n_frames), [0, 5])
        assert np.all(arr2 == arr3)

        # Trajectory
        traj = traj[:]
        arr = pt.distance(traj(stop=4), [0, 5])
        arr1 = pt.distance(traj(stop=4), [0, 5], n_frames=4)
        assert np.all(arr == arr1)

        arr2 = pt.distance(traj(stop=1000), [0, 5])
        arr3 = pt.distance(traj(stop=traj.n_frames), [0, 5])
        assert np.all(arr2 == arr3)

    def test_distance_with_dry_traj_and_PBC_topology(self):
        '''Situation: there is dry traj (no box) but Topology has box info
        '''
        traj = pt.iterload('data/dry_traj_with_PBC_top/strip.nc',
                           'data/dry_traj_with_PBC_top/strip.prmtop')
        assert traj.top.has_box(), 'Topology must have box for testing'

        correct_distance_with_image_True = pt.distance(traj,
                                                       ':8@OP2 :5@N1',
                                                       image=True)
        correct_distance_with_image_False = pt.distance(traj,
                                                        ':8@OP2 :5@N1',
                                                        image=False)
        state = pt.load_batch(traj, '''
        distance :8@OP2 :5@N1
        ''')
        state.run()
        expected_distance = [3.08030475, 2.68452183]

        aa_eq(correct_distance_with_image_False, expected_distance)
        aa_eq(correct_distance_with_image_True, [0., 0.])


class TestPairwiseDistance(unittest.TestCase):

    def test_pairwise(self):
        traj = pt.iterload('data/tz2.nc', 'data/tz2.parm7')
        distances = pt.pairwise_distance(traj, '@CA', '@CB')[0]

        ca_indices = pt.select_atoms('@CA', traj.top)
        cb_indices = pt.select_atoms('@CB', traj.top)
        known_shape = (traj.n_frames, len(ca_indices), len(cb_indices))
        assert known_shape == distances.shape, 'distance array shape'

        slow_distances = []
        for ca_i in ca_indices:
            for cb_i in cb_indices:
                slow_distances.append(pt.distance(traj, [ca_i, cb_i]))
        slow_distances = np.array(slow_distances).T
        aa_eq(slow_distances.flatten(), distances.flatten())


if __name__ == "__main__":
    unittest.main()
