#!/usr/bin/env python

from __future__ import print_function
import pytraj as pt


def test_0():
    traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")
