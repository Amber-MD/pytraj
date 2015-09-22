import unittest
from time import time
from pytraj.base import *
from pytraj import allactions
from pytraj import io as mdio
from pytraj import adict
from pytraj.datasets import cast_dataset


class TestActionList(unittest.TestCase):
    def test_run_0(self):
        # load traj
        farray = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        # creat datasetlist to hold distance data
        dsetlist = DatasetList()
        dflist = DataFileList()

        # creat ActionList to hold actions
        alist = ActionList()

        # creat TopologyList
        toplist = TopologyList()

        # add parm
        toplist.add_parm(farray.top)

        alist.add_action(allactions.Action_Distance(),
                         ArgList(":2@CA :3@CA out ./output/test_df.dat"),
                         toplist, dsetlist, dflist)
        alist.add_action(adict['dihedral'],
                         ":2@CA :3@CA :4@CA :5@CA out ./output/_dih.out",
                         toplist, dsetlist, dflist)

        #
        #print("test setup_actions")
        #print(alist.n_actions)

        # do checking
        alist.process(toplist[0])

        farray2 = Trajectory()
        frame0 = Frame()
        # testing how fast to do the actions

        t0 = time()
        # loop all frames
        # use iterator to make faster loop
        # don't use "for i in range(farray.size)"
        for frame in farray:
            # perform actions for each frame
            frame0 = frame.copy()
            alist.do_actions(frame0)
            # alist.do_actions(frame)

            # we need to keep the modified frame in farray2
            # farray2.append(frame)
            farray2.append(frame0)
        #print(time() - t0)


if __name__ == "__main__":
    unittest.main()
