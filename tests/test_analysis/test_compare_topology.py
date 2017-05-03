import pytraj as pt
from pytraj.analysis.topology_analysis import compare_topology
from pytraj.datafiles import load_tz2

from utils import fn

expected = """
#< tz2.parm7
#> tz2.parm7
# Atom types
# LJ params
# Bonds
# Bond Parameters
# Angles
# Angle Parameters
# Dihedrals
# Dihedral Parameters
"""


def test_compare_topology():
    # TODO: this is dummy test for coverage
    # wait until cpptraj has a proper test for this
    top = pt.load_topology(fn('tz2.parm7'))
    out = compare_topology(top, top)
    assert out.strip() == expected.strip()
