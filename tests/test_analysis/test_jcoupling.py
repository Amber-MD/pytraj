import os
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq, cpptraj_test_dir


def test_jcoupling():
    kfile = os.path.abspath(
        os.path.join(cpptraj_test_dir, "../dat/Karplus.txt"))
    traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))

    d1 = pt.jcoupling(traj, kfile=kfile)

    # load cpptraj
    txt = '''
    parm {}
    trajin {}
    jcoupling kfile {}
    '''.format(fn('tz2.parm7'), fn('tz2.nc'), kfile)
    cpptraj_out = pt.datafiles.load_cpptraj_output(txt)[1:]
    aa_eq(d1.values, cpptraj_out.values)
