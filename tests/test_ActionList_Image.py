import unittest
import pytraj as pt
from time import time
from pytraj.base import *
from pytraj import allactions
from pytraj import io as mdio
from pytraj.utils.check_and_assert import file_exist


class TestActionList(unittest.TestCase):
    def test_run_0(self):
        # load traj
        farray = Trajectory(
            filename="./data/tz2.truncoct.nc",
            top="./data/tz2.truncoct.parm7")[:2]
        fold = farray.copy()

        act = allactions.Action_Image()
        ptrajin = """                                                     
        center :2-11
        image center familiar com :6                                      
        """

        # create 'strip' action
        stripact = allactions.Action_Strip()

        # creat datasetlist to hold distance data
        dsetlist = DataSetList()
        dflist = DataFileList()

        # creat ActionList to hold actions
        alist = ActionList()

        # creat TopologyList
        toplist = TopologyList()

        # add parm
        toplist.add_parm(farray.top)

        # add two actions: Action_Strip and Action_Distance
        alist.add_action(
            allactions.Action_Center(), ArgList(":2-11"),
            top=toplist)
        alist.add_action(
            allactions.Action_Image(), ArgList("center familiar com :6"),
            top=toplist)

        #
        assert alist.n_actions == 2

        # do checking
        alist.process(toplist[0])

        farray2 = Trajectory()
        frame0 = Frame()
        # testing how fast to do the actions

        # loop all frames
        # use iterator to make faster loop
        # don't use "for i in range(farray.n_frames)"
        for frame in farray:
            # perform actions for each frame
            # we make a copy since we want to keep orginal Frame
            frame0 = frame.copy()
            alist.do_actions(frame0)
            # alist.do_actions(frame)

            # we need to keep the modified frame in farray2
            # farray2.append(frame)
            farray2.append(frame0)

        # make sure that Action_Strip does its job in stripping
        assert farray2.n_frames == farray.n_frames

        fsaved = mdio.iterload("./CpptrajTest/Test_Image/image4.crd.save",
                               "./data/tz2.truncoct.parm7")
        assert fsaved.n_frames == 2
        # make sure that pytraj reproduce cpptraj outputo
        # TODO : not yet. double-check


if __name__ == "__main__":
    unittest.main()
