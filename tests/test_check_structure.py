from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca

'''Aim: check segmentation fault for Topology that does not have bond info'''
class Test(unittest.TestCase):
    def test_0(self):
        traj = io.load_sample_data('tz2').to_mutable_trajectory()
        tbad = traj.copy()

        # make atom 0 and 1 very close
        tbad[0, 0] = 1.2
        tbad[1, 0] = 1.3

        check = adict['checkstructure']
        from pytraj import DataSetList
        dslist = DataSetList()
        check("", traj[0], top=traj.top, dslist=dslist)
        print (dslist)
        pyca.check_structure(tbad[0], top=tbad.top)


if __name__ == "__main__":
    unittest.main()
