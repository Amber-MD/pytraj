import numpy as np
import pytraj as pt
from pytraj.testing import aa_eq
from pytraj.testing import tempfolder
from pytraj.testing import cpptraj_test_dir

print("cpptraj_test_dir", cpptraj_test_dir)


def fn(name):
    return cpptraj_test_dir + '/Test_AtomMap/' + name


def get_content(fn):
    with open(fn) as fh:
        return fh.read()


def test_atom_map():
    start_fn = fn('start.mol2')
    xta_fn = fn('xtallig.mol2')
    saved_atommap_file = fn('atommap.dat.save')
    saved_xta_reordered_file = fn('reordered.pdb.save')
    saved_xta_reordered_traj = pt.iterload(saved_xta_reordered_file)
    saved_rmsd_file = fn('rmsd.dat.save')

    saved_data = get_content(saved_atommap_file)

    # using iterload to make sure trajectory coordinates won't be changed
    traj = pt.iterload(xta_fn)
    ref = pt.iterload(start_fn)

    # default
    data = pt.atom_map(traj, ref=ref)
    assert data[0] == saved_data

    # rmsfit
    with tempfolder():
        data = pt.atom_map(traj, ref=ref, rmsfit=True)
        mask_out = data[0]
        assert mask_out == saved_data
        saved_rmsd = np.loadtxt(saved_rmsd_file, skiprows=1).T[1]
        aa_eq(data[1], saved_rmsd)
