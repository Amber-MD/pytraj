import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.utils.check_and_assert import get_amber_saved_test_dir

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")
        f0 = traj[0].copy()
        f0saved = traj[0].copy()

        print (f0[:2])
        print (f0.has_box())
        act = adict['center']
        act(":1-13", f0, traj.top)

        print (f0[:2])

        saved_test = get_amber_saved_test_dir("Test_Center/centered.crd.save")
        if saved_test:
            print ("has saved test")
            saved_f = mdio.load(saved_test, "./data/tz2.truncoct.parm7")[0]
            print (saved_f[:2])

            # make sure we did the right thing
            assert_almost_equal(saved_f.coords, f0.coords)
            print ("OK")
        else:
            print ("can not find saved test")

if __name__ == "__main__":
    unittest.main()
