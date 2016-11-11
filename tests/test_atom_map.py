import numpy as np
import pytraj as pt
from pytraj.testing import tempfolder
from pytraj.testing import cpptraj_test_dir

def fn(name):
    return cpptraj_test_dir + '/Test_AtomMap/' + name 

def test_atom_map():
    start_fn = fn('start.mol2')
    xta_fn = fn('xtallig.mol2')
    saved_atommap_file = fn('atommap.dat.save')
    traj = pt.load(xta_fn)
    ref = pt.load(start_fn)

    with open(saved_atommap_file) as fh:
        saved_data = fh.read()

    data = pt.atom_map(traj, ref=ref)
    assert data == saved_data
