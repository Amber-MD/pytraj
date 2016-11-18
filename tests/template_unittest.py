#!/usr/bin/env python

from __future__ import print_function
import pytraj as pt
from utils import fn


def test_0():
    traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))
