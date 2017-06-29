#!/usr/bin/env python

from __future__ import print_function
import os
import pytraj as pt
from pytraj.testing import cpptraj_test_dir, aa_eq
from utils import fn


def test_fiximagebonds_default():
    pdb_fn = os.path.join(cpptraj_test_dir, 'Test_FixImagedBonds', 'MET.pdb')
    expected_rst = os.path.join(cpptraj_test_dir, 'Test_FixImagedBonds',
                                'fixed.rst7.save')

    traj = pt.load(pdb_fn)
    expected_traj = pt.load(expected_rst, top=traj.top)

    pt.fiximagedbonds(traj)
    aa_eq(expected_traj.xyz, traj.xyz)


def test_fiximagebonds_with_image():
    parm7 = os.path.join(cpptraj_test_dir, 'tz2.truncoct.parm7')
    traj_fn = os.path.join(cpptraj_test_dir, 'tz2.truncoct.nc')
    expected_crd = os.path.join(cpptraj_test_dir, 'Test_FixImagedBonds',
                                'unimage.crd.save')

    traj = pt.load(traj_fn, top=parm7)

    pt.image(traj, 'byatom')
    pt.fiximagedbonds(traj, ':1-13')
    traj.strip(':WAT')
    expected_traj = pt.load(expected_crd, top=traj.top)

    aa_eq(expected_traj.xyz, traj.xyz, decimal=3)
