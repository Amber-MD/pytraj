#!/usr/bin/env python

import pytraj as pt
from pytraj.testing import Timer
import mdtraj as md

fname, tname = 'GAAC3.1000frames.nc', 'GAAC.parm7'

@Timer()
def load_netcdf_(fname=fname, tname=tname):
    pt.io._load_netcdf(fname, tname)

@Timer()
def load_normal(fname=fname, tname=tname):
    pt.load(fname, tname)

@Timer()
def load_mdtraj(fname=fname, tname=tname):
    md.load_netcdf(fname, top=tname)

for _ in range(5):
    print(_)
    load_netcdf_()
    load_normal()
    load_mdtraj()
