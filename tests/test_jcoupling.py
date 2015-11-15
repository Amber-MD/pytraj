import os
import unittest
import pytraj as pt
from pytraj.testing import aa_eq, cpptraj_test_dir


class Test(unittest.TestCase):

    def test_0(self):
        kfile = os.path.abspath(os.path.join(cpptraj_test_dir,
                                             "Test_Jcoupling", "Karplus.txt"))
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")

        d1 = pt.jcoupling(traj, kfile=kfile)

        # load cpptraj
        txt = '''
        parm ./data/tz2.parm7
        trajin data/tz2.nc
        jcoupling kfile %s
        ''' % kfile
        cpptraj_out = pt.datafiles.load_cpptraj_output(txt)[1:]
        aa_eq(d1.values, cpptraj_out.values)


if __name__ == "__main__":
    unittest.main()
