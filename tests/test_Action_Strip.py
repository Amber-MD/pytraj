import unittest
from pytraj.base import *
from pytraj.actions.Action_Strip import Action_Strip
from pytraj import allactions
from pytraj.actions import Action
from pytraj.TrajectoryIterator import TrajectoryIterator
from pytraj.decorators import no_test

print(dir(Action_Strip()))

farray = Trajectory(top=Topology("./data/Tc5b.top"), filename='data/md1_prod.Tc5b.x')

class TestStrip(unittest.TestCase):
    #@no_test
    def test_master(self):
        top = Topology("./data/Tc5b.top")
        newtop = top.copy()
        frame0 = farray[0].copy()
        newframe = Frame()
        act = Action_Strip()
        act_surf = allactions.Action_Surf()

        for i in range(1):
            current_frame = farray[i]
            newframe = Frame()
            dslist = DataSetList()
            act(command="strip !@CA", 
                       top=top, 
                       dslist=dslist,
                       current_frame=frame0, 
                       new_frame=newframe, 
                       new_top=newtop)

            act_surf(command="@CA", 
                       top=top, 
                       dslist=dslist,
                       current_frame=farray)

            act_surf(command="@H=", 
                       top=top, 
                       dslist=dslist,
                       current_frame=farray)

        print(newtop)
        print(newframe)
        print(dslist.is_empty())
        print(dslist.size)
        print(dir(dslist[0]))
        dcast = cast_dataset(dslist[0], dtype='general')
        print(dcast)
        print(dcast.size)
        print(dcast[:20])

    #@no_test
    def test_0(self):
        print("newtop")
        farray0 = farray.copy()
        newtop = farray0.top.copy()
        oldtop = farray0.top
        
        toplist = TopologyList()
        toplist.add_parm(newtop)
        dslist = DataSetList()
        dflist = DataFileList()
        
        stripact = Action_Strip()
        #stripact.read_input("strip !@CA", toplist)
        stripact.read_input("strip !@CA", oldtop)
        stripact.process(oldtop, newtop)

        frame0 = Frame(farray.top.n_atoms)
        stripact.do_action(farray[0], frame0)
        print(frame0.size)
        assert frame0.size == 60
        assert frame0.size != farray0[0].size

        stripact.do_action(farray[0], farray[0])
        print("after stripping", farray[0].size)

        frame0_view = farray0[0]
        stripact.do_action(farray0[0], frame0_view)
        print("after stripping", farray0[0].size)
        print(frame0_view.n_atoms)
        print(newtop.n_atoms)

    def test_1(self):
        from pytraj import adict
        farray0 = farray.copy()
        act = adict['strip']
        tmpframe = Frame()
        newf = Trajectory()
        newf.top = farray0.top.strip_atoms("!@CA", copy=True)

        act("!@CA", farray0, farray0.top.copy(), new_frame=tmpframe)
        print (farray0[0].n_atoms)
        print (tmpframe.n_atoms)

if __name__ == "__main__":
    unittest.main()
