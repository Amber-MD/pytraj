import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj import adict

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0 = traj[0]
        f0cp = f0.copy()
        act = adict['strip']
        act.read_input('!@CA', traj.top.copy())
        act.process(traj.top.copy())
        act.do_action(f0)
        newf = Frame()
        act.do_action(f0, newf)
        print (f0)
        print (newf)
        print (newf[1])
        act.do_action(traj, newf)
        flast = traj[-1]
        flast.strip_atoms("!@CA", traj.top)
        print (newf)
        print (newf[1])
        print (flast[1])

        newf = Frame()
        act2 = adict['strip']
        act2("!@CA", f0cp, traj.top, new_frame=newf)
        # FIXME : newf must has 20 atoms, but it has 284 atoms
        print (newf)
        print (newf[1])

    def test_1(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0 = traj[0]
        f0cp = f0.copy()
        act = adict['strip']
        act.read_input('!@CA', traj.top.copy())
        act.process(traj.top.copy())

        farray = FrameArray()
        newf = Frame()
        for i, frame in enumerate(traj):
            act.do_action(frame, newf)
            farray.append(newf)

        for frame in farray:
            assert frame.n_atoms == 20

if __name__ == "__main__":
    unittest.main()
