from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    def test_0(self):
        from glob import glob
        flist = glob("data/Test_RemdTraj/rem.nc.0*")
        tlist = mdio._iterload_from_filelist(flist, 
                "./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7",
                force_load=True)
        for t in tlist:
            print (t.filename)

if __name__ == "__main__":
    unittest.main()
