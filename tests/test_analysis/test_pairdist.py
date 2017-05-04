#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.testing import aa_eq, tempfolder

# local
from utils import fn


def test_pairdist():
    traj = pt.iterload(fn("tz2.crd"), fn("tz2.parm7"))

    for (mask, delta) in [('*', 0.1), ('@CA', '0.2')]:
        data = pt.pairdist(traj, delta=delta, mask=mask, dtype='ndarray')
        data0 = data[0].T
        data1 = data[1].T

        txt = '''
        parm {0}
        trajin {1}
        pairdist out test.out mask {2} delta {3}
        '''.format(fn('tz2.parm7'), fn('tz2.crd'), mask, delta)

        # take data, skip DatasetTopology
        state = pt.load_cpptraj_state(txt)
        with tempfolder():
            state.run()
        cpp_data = state.data[1:].values
        cpp_distance, cpp_Pr = cpp_data[0].T
        _, cpp_std = cpp_data[1].T

        aa_eq(data0[0], cpp_distance)  # distance
        aa_eq(data1[0], cpp_distance)  # distance
        aa_eq(data0[1], cpp_Pr, decimal=2)  # Pr
        aa_eq(data1[1], cpp_std, decimal=2)  # std
