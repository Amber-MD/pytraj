#!/usr/bin/env python

import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
from pytraj.utils.get_common_objects import super_dispatch


class TestSuperDispatch(unittest.TestCase):
    def setUp(self):
        self.traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))

    def test_naive(self):
        # make sure to convert int to Frame
        def func_convert_ref(traj, ref=0, top=None):
            assert isinstance(ref, pt.Frame)

        func = super_dispatch()(func_convert_ref)
        func(self.traj, ref=-2)

        # make sure to insert correct Topology
        def func_convert_top(traj, top=None):
            assert isinstance(top, pt.Topology)

        func = super_dispatch()(func_convert_top)
        func(self.traj, top=None)

        # make sure to convert array to Amber mask
        def func_convert_mask_array(traj, top=None, mask=None):
            assert isinstance(mask, str)

        func = super_dispatch()(func_convert_mask_array)
        func(self.traj, mask=[0, 3, 7])

        # test all: top, mask, ref
        def func_all_3(traj, mask='', ref=0, top=None):
            assert isinstance(mask, str)
            assert isinstance(ref, pt.Frame)
            assert isinstance(top, pt.Topology)

        func = super_dispatch()(func_all_3)
        func(self.traj, ref=3)

        super_dispatch()(func_all_3)(self.traj, ref=3)
        # specify nothing
        super_dispatch()(func_all_3)(self.traj)

    def testsuper_dispatch(self):
        traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))

        funclist = [pt.radgyr, pt.molsurf]
        for func in funclist:
            mask = '@CA'
            atom_indices = pt.select_atoms(mask, traj.top)
            # mask
            aa_eq(func(traj, mask=mask), func(traj, mask=atom_indices))
            # specify traj=traj
            aa_eq(func(traj=traj, mask=mask), func(traj, mask=atom_indices))

            # frame_indices with mask
            frame_indices = [0, 5, 8]
            aa_eq(
                func(traj[frame_indices], mask=mask),
                func(traj, mask=atom_indices, frame_indices=frame_indices))


if __name__ == "__main__":
    unittest.main()
