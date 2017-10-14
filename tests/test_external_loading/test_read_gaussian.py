from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq

try:
    import cclib
    has_cclib = True
except ImportError:
    has_cclib = False


@unittest.skipUnless(has_cclib, 'skip: does not have cclib')
def test_read_gaussian_output():
    filename = "./data/gaussian/GF2.log"
    gau = cclib.parser.Gaussian(filename)
    go = gau.parse()

    traj = pt.tools.read_gaussian_output(filename, "./data/gaussian/GF2.pdb")
    aa_eq(traj.xyz, go.atomcoords)
