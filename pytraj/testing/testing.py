from __future__ import absolute_import, print_function
import os

from pytraj.datafiles.load_samples import load_sample_data
from pytraj.utils import eq
from pytraj.utils import Timer
from pytraj.utils import tempfolder
from pytraj.utils import duplicate_traj
from numpy.testing import assert_almost_equal
from functools import partial
import numpy as np

aa_eq = partial(assert_almost_equal, decimal=4)

__all__ = [
    'load_sample_data', 'eq', 'aa_eq', 'duplicate_traj', 'Timer', 'tempfolder',
    'amberhome', 'cpptraj_test_dir', 'get_fn', 'get_remd_fn',
    'assert_equal_topology', 'assert_equal_dict', 'load_cpptraj_reference_data'
]

# find cpptraj test dir
# find in CPPTRAJHOME first
# if not exist CPPTRAJHOME, find in AMBERHOME
# if not exist CPPTRAJHOME and AMBERHOME, find in ../cpptraj/ folder
# (supposed we run pytraj's test in pytraj/tests folder)

# if you are adding more test files to cpptraj master branch on github,
# you should export CPPTRAJHOME

cpptrajhome = os.environ.get('CPPTRAJHOME', '')
amberhome = os.environ.get('AMBERHOME', '')

if cpptrajhome:
    cpptraj_test_dir = os.path.join(cpptrajhome, 'test')
else:
    if amberhome:
        cpptraj_test_dir = os.path.join(amberhome, 'AmberTools', 'src',
                                        'cpptraj', 'test')
    else:
        cpptrajhome = ''
        amberhome = ''
        cpptraj_test_dir = ''

DEFAULT_PATH = os.path.dirname(__file__) + "/../../cpptraj/test"

if os.path.exists(DEFAULT_PATH):
    cpptraj_test_dir = os.path.abspath(DEFAULT_PATH)

header_doc = '''
import pytraj as pt
import pytraj.all_actions as pyca
traj = pt.load_sample_data("tz2")
'''


def assert_equal_topology(top, new_top, traj=None):
    import pytraj as pt
    assert new_top.n_atoms == top.n_atoms, 'same n_atoms'
    assert new_top.n_residues == top.n_residues, 'same n_residues'
    assert new_top.n_mols == top.n_mols, 'same n_mols'
    # there are inverted bond indices [5292 5291] vs [5291 5292]
    # so use distance to assert
    if traj is not None:
        aa_eq(
            pt.distance(traj, new_top.bond_indices),
            pt.distance(traj, top.bond_indices))
        # same for dihedral_indices
        aa_eq(
            pt.dihedral(traj, new_top.dihedral_indices),
            pt.dihedral(traj, top.dihedral_indices))
    aa_eq(new_top.dihedral_indices, top.dihedral_indices)
    aa_eq(new_top.mass, top.mass)
    aa_eq(new_top.charge, top.charge)
    aa_eq(new_top.box.values, top.box.values)

    assert [res.name for res in top.residues] == [
        res.name for res in new_top.residues
    ], 'equal resnames'
    assert [atom.name for atom in top.atoms] == [
        atom.name for atom in new_top.atoms
    ], 'equal resnames'

    for res, res_new in zip(top.residues, new_top.residues):
        assert res.first_atom_index == res_new.first_atom_index, 'first atom'
        assert res.last_atom_index == res_new.last_atom_index, 'last atom'

def assert_equal_dict(dict1, dict2):
    """Assert two dict are equal (for testing purpose)
    Assume the values are numpy arrays
    """
    assert dict1.keys() == dict2.keys()
    for key in dict1.keys():
        val1 = dict1[key]
        val2 = dict2[key]
        if val1.dtype == np.dtype("<U9") or val1.dtype == np.dtype("object"):
            assert len(val1) == len(val2)
            assert all(v1 == v2 for v1, v2 in zip(val1, val2))
        else:
            aa_eq(val1, val2)


def get_fn(txt):
    '''get absolute path for trajectory and topology samples. Legit text = 'ala3', 'tz2',
    'rna'. Mostly for testing purpose.

    Examples
    --------
    >>> # get trajectory file name
    >>> fname = get_fn('tz2')[0]
    >>> fname.split('/')[-1]
    'tz2.ortho.nc'
    '''
    from pytraj import load_sample_data
    traj = load_sample_data(txt)
    return traj.filename, traj.top.filename


def get_remd_fn(txt):
    """

    >>> fnlist, tn = get_remd_fn('remd_ala2')
    """
    from pytraj import load_sample_data
    traj = load_sample_data(txt)
    return traj.filelist, traj.top.filename


def load_cpptraj_reference_data(test_name, data_file, column=1):
    """Load reference data from cpptraj test files

    Parameters
    ----------
    test_name : str
        Name of the cpptraj test directory (e.g., 'Test_Distance', 'Test_RMSD')
    data_file : str
        Name of the reference data file (e.g., 'dist.dat.save', 'rmsd.dat.save')
    column : int, optional
        Column index to extract (0-based), default is 1 (second column)

    Returns
    -------
    np.ndarray or None
        Array of reference values, or None if file not found

    Examples
    --------
    >>> # Load distance reference data
    >>> distances = load_cpptraj_reference_data('Test_Distance', 'dist.dat.save')
    >>> # Load RMS fluctuations (5th column)
    >>> fluct = load_cpptraj_reference_data('Test_Analyze_Modes', 'fluct.dat.save', column=4)
    """
    if not cpptraj_test_dir:
        return None

    cpptraj_ref_file = os.path.join(cpptraj_test_dir, test_name, data_file)

    if not os.path.exists(cpptraj_ref_file):
        return None

    values = []
    with open(cpptraj_ref_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) > column:
                try:
                    values.append(float(parts[column]))
                except ValueError:
                    continue

    return np.array(values) if values else None


if __name__ == "__main__":
    print(amberhome)
    print(cpptraj_test_dir)
