from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = traj.search_hbonds()
        d0 = dslist[0]
        assert d0 in dslist
        d0cp = d0.copy()
        assert d0cp not in dslist

        # FIXME: assert failed
        dslist2 = dslist.__class__()
        # append with copying
        dslist2.append(d0cp)
        print (dslist2)
        assert d0cp not in dslist2
        dslist3 = dslist.__class__()
        dslist3.append(d0cp, copy=False)
        assert d0cp in dslist3

if __name__ == "__main__":
    unittest.main()
