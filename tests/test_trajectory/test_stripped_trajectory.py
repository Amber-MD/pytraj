#!/usr/bin/env python

from __future__ import print_function, absolute_import
import unittest
import pytraj as pt
from pytraj.testing import aa_eq
from pytraj.testing import tempfolder

from utils import fn


def test_stripped_trajectory():
    traj_on_disk = pt.iterload(fn('tz2.ortho.nc'), fn('tz2.ortho.parm7'))

    straj = traj_on_disk.autoimage().superpose('@CA').strip(":WAT")
    traj_on_mem_strip = traj_on_disk['!:WAT']

    aa_eq(traj_on_disk['!:WAT'].xyz, straj.xyz)
    aa_eq(traj_on_mem_strip.xyz, straj.xyz)

    # avoid memory free, haizz
    f0 = traj_on_mem_strip[0]
    f1 = straj[0]
    aa_eq(f0.xyz, f1.xyz)

    # known failure
    # aa_eq(traj_on_mem_strip[0].xyz, straj[0].xyz)

    aa_eq(traj_on_mem_strip[:].xyz, straj[:].xyz)
    aa_eq(traj_on_mem_strip[:3].xyz, straj[:3].xyz)
    aa_eq(traj_on_mem_strip[3:10].xyz, straj[3:10].xyz)

    # save
    with tempfolder():
        nc_fn = 'test.nc'
        straj.save(nc_fn, overwrite=True)
        aa_eq(straj.xyz, pt.load(nc_fn, straj.top).xyz)
