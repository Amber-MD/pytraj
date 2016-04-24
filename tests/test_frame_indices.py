#!/usr/bin/env python

from __future__ import print_function
import unittest
from pytraj import *
import numpy as np
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.compat import string_types
from pytraj.hbond_analysis import DatasetHBond


class TestFrameIndices(unittest.TestCase):

    def setUp(self):
        self.traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
        self.traj_ortho = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        self.traj_nu = pt.iterload('data/Test_NAstruct/adh026.3.pdb')

    def test_frame_indices_from_yield(self):
        '''extensive and seperated testsing
        '''
        traj = self.traj

        def gen_int():
            for i in range(0, 10, 2):
                yield i

        for idx, frame in enumerate(pt.iterframe(traj,
                                                 frame_indices=gen_int())):
            pass

        assert idx == len(range(0, 10, 2)) - 1
        aa_eq(frame.xyz, traj[8].xyz)

    def test_frame_indices_for_function(self):
        traj = self.traj

        pdict = pt.__dict__
        funclist = list(set(pdict[key]
                            for key in dir(pt)
                            if hasattr(pdict[key], '_issuper_dispatched')))

        frame_indices = [0, 2]

        # remove 'calc_jcoupling' since does not have kfile on travis
        # remove energy_decomposition since does not have sander
        # remove center, why?
        # remove search_neighbors, why? (got messup with Frame memory owner)
        excluded_fn = [jcoupling, volmap,
                       center, search_neighbors,
                       atomiccorr, autoimage, closest,
                       volume, superpose, randomize_ions,
                       check_structure,
                       align_principal_axis, ]
        func_nu = [
            calc_epsilon, calc_alpha, calc_zeta, calc_beta, calc_nu1, calc_nu2,
            calc_delta, calc_chin,
            calc_gamma, ]

        # default mask, default ref
        for func in funclist:
            if func not in excluded_fn:
                if func is pt.calc_multivector:
                    data_0 = func(traj,
                                  resrange='1-6',
                                  names='C N',
                                  frame_indices=frame_indices)
                    data_1 = func(traj[frame_indices],
                                  resrange='1-6',
                                  names='C N')
                elif func is pt.volmap:
                    # use water
                    data_0 = func(self.traj_ortho,
                                  mask=':WAT@O',
                                  grid_spacing=(0.2, 0.2, 0.2),
                                  centermask='!:1-13',
                                  frame_indices=frame_indices)
                    data_1 = func(self.traj_ortho[frame_indices],
                                  mask=':WAT@O',
                                  centermask='!:1-13',
                                  grid_spacing=(0.2, 0.2, 0.2))
                elif func in func_nu:
                    data_0 = func(self.traj_nu, frame_indices=frame_indices)
                    data_1 = func(self.traj_nu[frame_indices])
                else:
                    data_0 = func(traj, frame_indices=frame_indices)
                    data_1 = func(traj[frame_indices])

                if isinstance(data_0, np.ndarray):
                    aa_eq(data_0, data_1)
                elif isinstance(data_0, pt.DatasetList):
                    for arr0, arr1 in zip(data_0, data_1):
                        # do each element in case we can not convert DatasetList to
                        # ndarray
                        if not isinstance(arr0[0], string_types):
                            aa_eq(arr0.values, arr1.values)
                elif isinstance(data_0, DatasetHBond):
                    aa_eq(data_0.data.values, data_1.data.values)
                elif isinstance(data_0, (list, tuple)):
                    # dssp
                    aa_eq(data_0[-1].values, data_1[-1].values)
                else:
                    raise RuntimeError(
                        'must return ndarray or DatasetList or DatasetHBond')

        # test excluded fns
        aa_eq(pt.atomiccorr(traj[frame_indices], '@CA'),
              pt.atomiccorr(traj, '@CA', frame_indices=frame_indices))

        # align_principal_axis
        indices = [0, 3, 7]
        t0 = self.traj[:]
        t1 = self.traj[indices]
        pt.align_principal_axis(t0, frame_indices=indices)
        pt.align_principal_axis(t1)
        aa_eq(t0[indices].xyz, t1.xyz)
        # make sure that other frames are not affected
        other_indices = list(set(range(10)) - set(indices))
        aa_eq(self.traj[other_indices].xyz, t0[other_indices].xyz)


if __name__ == "__main__":
    unittest.main()
