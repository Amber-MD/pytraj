import unittest
from pytraj.common_actions import calc_molsurf
from pytraj.base import *
from pytraj.actions.Action_Surf import Action_Surf
from pytraj.actions import Action
from pytraj.TrajReadOnly import TrajReadOnly
from pytraj import cast_dataset
from pytraj import adict

print(dir(Action_Surf()))

farray = FrameArray(top=Topology("./data/Tc5b.top"), filename='data/md1_prod.Tc5b.x')

class TestSurf(unittest.TestCase):
    def test_0(self):
        print("newtop")

        farray0 = farray.copy()
        newtop = farray0.top.copy()
        oldtop = farray0.top

        d0 = calc_molsurf("@CA", farray[0], farray.top)
        print (d0[:])
        adict['molsurf'].help()
        
        toplist = TopologyList()
        toplist.add_parm(newtop)
        dslist = DataSetList()
        dflist = DataFileList()
        
        act = Action_Surf()
        act.read_input("@CA", oldtop, dslist=dslist)
        act.process(oldtop, newtop)
        frame0 = Frame(farray.top.n_atoms)
        act.do_action(farray[0], frame0)
        print (dslist[0][:])
        act.help()

if __name__ == "__main__":
    unittest.main()
