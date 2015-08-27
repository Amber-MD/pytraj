import unittest
import pytraj as pt
from time import time
from pytraj import allactions
from pytraj import adict
from pytraj.base import *


class TestActionList(unittest.TestCase):
    def test_run_0(self):
        # load traj
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = DataSetList()
        dflist = DataFileList()

        # creat ActionList to hold actions
        alist = ActionList()
        # add two actions: Action_Dihedral and Action_Distance
        alist.add_action(adict['distance'],
                         ":2@CA :10@CA out ./output/_dist.out", traj.top,
                         dslist, dflist)
        alist.add_action(adict['dihedral'],
                         ":2@CA :3@CA :4@CA :5@CA out ./output/_dih.out",
                         traj.top, dslist, dflist)

        # using string for action 'dssp'
        alist.add_action('dssp', "out ./output/_dssp_alist.out", traj.top,
                         dslist, dflist)
        alist.add_action('matrix', "out ./output/_mat_alist.out", traj.top,
                         dslist, dflist)
        # does not work with `strip` (output traj have the same n_atoms as originl traj)
        #alist.add_action("strip", "!CA", traj.top)
        alist.add_action("outtraj", "./output/test_trajout.nc", traj.top)
        alist.do_actions([traj[[0, 1]], traj, traj.chunk_iter(chunksize=4,
                                                              stop=8),
                          traj.iterframe()])
        Nframes = 1 + 1 + traj.n_frames + 8 + traj.n_frames
        dflist.write_all_datafiles()
        traj2 = pt.iterload("./output/test_trajout.nc", traj.top)
        assert traj2.n_frames == Nframes

    def test_run_1(self):
        # load traj
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = DataSetList()
        dflist = DataFileList()

        # creat ActionList to hold actions
        alist = ActionList()
        alist.add_action(adict['distance'],
                         ":2@CA :10@CA out ./output/_dist.out", traj.top,
                         dslist, dflist)
        alist.do_actions(traj.iterframe())
        assert dslist.size == 1
        assert dslist[0].size == traj.size

    def test_run_1(self):
        # load traj
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = DataSetList()
        dflist = DataFileList()

        # creat ActionList to hold actions
        alist = ActionList()
        alist.add_action(adict['distance'],
                         ":2@CA :10@CA out ./output/_dist.out", traj.top,
                         dslist, dflist)
        alist.do_actions([traj.chunk_iter()])
        assert dslist.size == 1
        assert dslist[0].size == traj.size


if __name__ == "__main__":
    unittest.main()
