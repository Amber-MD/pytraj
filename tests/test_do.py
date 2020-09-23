#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.testing import aa_eq

# local
from utils import fn


def test_pytraj_do():
    traj = pt.iterload(fn("tz2.ortho.nc"), fn("tz2.ortho.parm7"))

    ref0 = pt.autoimage(traj[0], top=traj.top)
    ref1 = pt.autoimage(traj[1], top=traj.top)

    data = pt.tools.dict_to_ndarray(
        pt.compute(
            [
                'autoimage', 'radgyr nomax @CA', 'rms refindex 0',
                'rms refindex 1'
            ],
            traj,
            ref=[ref0, ref1]))

    t0 = pt.autoimage(traj[:])
    aa_eq(pt.radgyr(t0, '@CA'), data[0])
    aa_eq(pt.rmsd(t0, ref=ref0), data[1])
    aa_eq(pt.rmsd(t0, ref=ref1), data[2])
    aa_eq(pt.compute(pt.radgyr, t0, '@CA'), data[0])
