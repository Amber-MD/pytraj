#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestTransformation(unittest.TestCase):

    def test_transforming_trajectory(self):
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        t0 = traj[:]

        t0.transform(['autoimage', 'center origin', 'rotate x 30',
                      'scale x 1.2', 'translate x 1.3'])
        t1 = traj[:].autoimage().center('origin').rotate('x 30').scale('x 1.2').translate('x 1.3')
        aa_eq(t0.xyz, t1.xyz)

        t2 = traj[:]
        pt.transform(t2, ['autoimage', 'center origin', 'rotate x 30',
                          'scale x 1.2', 'translate x 1.3'])
        aa_eq(t1.xyz, t2.xyz)

        # another way
        t3 = traj[:]
        t3.autoimage().center('origin')
        pt.rotate(t3, 'x 30')
        pt.scale(t3, 'x 1.2')
        pt.translate(t3, 'x 1.3')
        aa_eq(t3.xyz, t2.xyz)

    def test_align_principal_axis(self):
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        t0 = traj[:]
        t1 = traj[:]
        aa_eq(t0.align_principal_axis().xyz, pt.align_principal_axis(t1).xyz)

        # compare to cpptraj
        cm = '''
        principal dorotation
        createcrd mycrd
        '''
        state = pt.load_cpptraj_state(cm, traj)
        state.run()
        aa_eq(state.data['mycrd'].xyz, t1.xyz)

if __name__ == "__main__":
    unittest.main()
