import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.actions.Action_Energy import Action_Energy

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = Action_Energy()
        print (act)
        dslist = DataSetList()
        act.read_input("bond angle dihedral", traj.top, dslist=dslist)
        act.help()
        act.process(traj.top)
        act.do_action(traj)
        print (dslist.size)
        size = dslist.size
        for i in range(size):
            print (dslist[i])
            print (dslist[i].name)
        print (dslist[0][:])
        d0 = dslist[0]

if __name__ == "__main__":
    unittest.main()
