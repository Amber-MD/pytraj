import unittest
import pytraj as pt
from pytraj.base import *
from pytraj.actions.CpptrajActions import Action_Strip
from pytraj import allactions
from pytraj.actions import Action
from pytraj.TrajectoryIterator import TrajectoryIterator
from pytraj.decorators import no_test


farray = Trajectory(
    top=Topology("./data/Tc5b.top"),
    filename='data/md1_prod.Tc5b.x')


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

        dcast = cast_dataset(dslist[0], dtype='general')

    def test_0(self):
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
        assert frame0.size == 60
        assert frame0.size != farray0[0].size

        stripact.do_action(farray[0], farray[0])

        frame0_view = farray0[0]
        stripact.do_action(farray0[0], frame0_view)

    def test_1(self):
        from pytraj import adict
        farray0 = farray.copy()
        act = adict['strip']
        tmpframe = Frame()
        newf = Trajectory()
        newf.top = farray0.top.strip_atoms("!@CA", copy=True)

        act("!@CA", farray0, farray0.top.copy(), new_frame=tmpframe)


if __name__ == "__main__":
    unittest.main()
