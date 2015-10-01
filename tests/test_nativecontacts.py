from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        txt = pt.datafiles.tc5b_trajin + '''
        nativecontacts
        '''

        dslist = pt.native_contacts(traj, top=traj.top)
        dslist1 = pt.native_contacts(traj, top=traj.top, ref=-1)
        cpp = pt.datafiles.load_cpptraj_output(txt, dtype='state')
        # remove DatasetTopology
        cpp.data.remove_set(cpp.data[0])
        aa_eq(dslist.values, cpp.data.values)

if __name__ == "__main__":
    unittest.main()
