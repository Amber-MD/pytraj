import unittest
import pytraj as pt
from pytraj.testing import aa_eq
from pytraj.compat import zip
from pytraj import adict


class TestAutoImage(unittest.TestCase):

    def test_1(self):
        traj = pt.iterload("./data/tz2.truncoct.nc",
                           "./data/tz2.truncoct.parm7")
        f0 = traj[0]
        f0cp = f0.copy()
        adict['autoimage']("", f0, traj.top)
        fsaved = pt.iterload("./data/tz2.truncoct.autoiamge.save.r",
                             "./data/tz2.truncoct.parm7")[0]
        aa_eq(fsaved.xyz, f0.xyz, decimal=3)

    def test_2(self):
        from pytraj.all_actions import do_autoimage
        # test do_autoimage
        traj = pt.iterload("./data/tz2.truncoct.nc",
                           "./data/tz2.truncoct.parm7")
        f0 = traj[0]
        f0cp = f0.copy()
        do_autoimage(traj=f0, top=traj.top)

        fsaved = pt.iterload("./data/tz2.truncoct.autoiamge.save.r",
                             "./data/tz2.truncoct.parm7")[0]
        aa_eq(fsaved.xyz, f0.xyz, decimal=3)


class TestGeometry(unittest.TestCase):

    def test_radgyr(self):
        traj = pt.iterload(top="./data/Tc5b.top",
                           filename='data/Tc5b.x', )
        txt = '''
        parm ./data/Tc5b.top
        trajin ./data/Tc5b.x
        radgyr @CA nomax
        radgyr nomax
        radgyr !@H= nomax
        '''

        # exclude DatasetTopology
        data = pt.datafiles.load_cpptraj_output(txt)[1:]
        for mask, out in zip(['@CA', '', '!@H='], data):
            aa_eq(pt.radgyr(traj, mask), out)


if __name__ == "__main__":
    unittest.main()
