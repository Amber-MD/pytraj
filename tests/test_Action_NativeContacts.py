from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        # TODO : assert
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        txt = pt.datafiles.tc5b_trajin + '''
        nativecontacts
        '''
        dslist = pt.native_contacts(traj, top=traj.top)
        cpp = pt.datafiles.load_cpptraj_output(txt)

        aa_eq(dslist.values, cpp.values)


if __name__ == "__main__":
    unittest.main()
