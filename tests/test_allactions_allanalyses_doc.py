import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj import info, adict
from pytraj import allactions
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_action(self):
        info(adict["pairdist"])
        print ()
        info(adict["dihedral"])

        dslist = DataSetList()

        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = allactions.Action_Dihedral()
        act.read_input("dihedral :4@C :5@N :5@CA :5@C range360  mass ", 
                       current_top=traj.top, dslist=dslist)
        print(traj.top)
        act.process(current_top=traj.top)
        print(act.process.__doc__)

        for i, frame in enumerate(traj):
            act.do_action(frame)

        d0 = cast_dataset(dslist[0])
        print(d0[:10])

    def test_action_2(self):
        pass

if __name__ == "__main__":
    unittest.main()
