#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
import numpy as np
from pytraj.decorators import _register_pmap
from pytraj.utils import eq, aa_eq


class TestFun(unittest.TestCase):

    def test_fun(self):
        '''for anything that does not need serious testing
        '''
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")

        state = pt.load_pipeline(traj, ['distance :3 :2', ]).compute()


@_register_pmap
def method(traj):
    # make global function (to be pickable)
    x = []
    for frame in traj:
        frame.xyz *= 2
        x.append(np.array(frame.xyz[0]))
    return x


class TestUserFunction(unittest.TestCase):

    def test_user_defined_method(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
        for n_cores in [3, 5, 7]:
            data = pt.pmap(method, traj, n_cores=n_cores)
            joint_data = pt.tools.flatten(x[1] for x in data)
        aa_eq(joint_data, pt.tools.flatten(method(traj)))


if __name__ == "__main__":
    unittest.main()
