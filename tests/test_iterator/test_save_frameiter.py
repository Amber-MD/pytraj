from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq


def test_save_frameiter(tmpdir):
    with tmpdir.as_cwd():
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))

        traj(0, 8, 2, mask='@CA').save('dummy_test0.nc', overwrite=True)
        pt.write_traj(
            'dummy_test1.nc', traj(0, 8, 2, mask='@CA'), overwrite=True)

        new_top = traj.top._get_new_from_mask('@CA')
        t0 = pt.iterload('dummy_test0.nc', new_top)
        t1 = pt.iterload('dummy_test1.nc', new_top)
        aa_eq(t0.xyz, t1.xyz)
