#!/usr/bin/env python
from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
from mock import patch


def test_nativecontacts(tmpdir):
    with tmpdir.as_cwd():
        traj_fname = fn('DPDP.nc')
        top_fname = fn('DPDP.parm7')
        traj = pt.iterload(traj_fname, top_fname)

        c_raw = pt.compute("""
        nativecontacts name NC3 :1-21&!@H= byresidue out nc.all.res.dat distance 3.0
        """, traj)

        data = pt.nativecontacts(
            traj,
            mask=':1-21&!@H=',
            byres=True,
            distance=3.0,
            options='name NC3',
            dtype='dict')
        print(c_raw, data)
        aa_eq(c_raw['NC3[native]'], data['NC3[native]'])
