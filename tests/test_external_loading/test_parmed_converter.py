#!/usr/bin/env python

import pytest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
from pytraj.sandbox import to_parmed

pytestmark = pytest.mark.skip(reason="Skipping the whole file for now")

try:
    import parmed as pmd
except ImportError:
    pmd = None

@pytest.fixture
def traj():
    return pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))

@pytest.fixture
def parm(traj):
    return pmd.load_file(traj.top.filename) if pmd else None

def parm_eq(parm, parm2):
    for atom, atom2 in zip(parm.atoms, parm2.atoms):
        assert atom.name == atom2.name, 'equal name'
        assert atom.type == atom2.type, 'equal type'
        assert atom.mass == atom2.mass, 'equal mass'
        assert atom.atomic_number == atom2.atomic_number, 'equal atomic_number'
        assert atom.residue.name == atom2.residue.name, 'residue name'
        aa_eq(atom.charge, atom2.charge, decimal=4)

@pytest.mark.skipif(pmd is None, reason="Must install ParmEd")
def test_parmed_converter(traj, parm):
    parm_eq(parm, to_parmed(traj))
    parm_eq(parm, traj.top.to_parmed())