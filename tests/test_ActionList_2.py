
import unittest
from time import time
from pytraj.base import *
from pytraj import allactions
from pytraj import io as mdio
from pytraj import adict

class TestActionList(unittest.TestCase):
    def test_run_0(self):
        # load traj
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        print (traj)
        dslist = DataSetList()
        dflist = DataFileList()
    
        # creat ActionList to hold actions
        alist = ActionList()
        # add two actions: Action_Dihedral and Action_Distance
        alist.add_action(adict['distance'], ":2@CA :10@CA out ./output/_dist.out", 
                         traj.top, dslist, dflist)
        alist.add_action(adict['dihedral'], ":2@CA :3@CA :4@CA :5@CA out ./output/_dih.out", 
                         traj.top, dslist, dflist)

        # using string for action 'dssp'
        alist.add_action('dssp', "out ./output/_dssp_alist.out", 
                         traj.top, dslist, dflist)
        alist.add_action('matrix', "out ./output/_mat_alist.out", 
                         traj.top, dslist, dflist)
        # does not work with `strip` (output traj have the same n_atoms as originl traj)
        #alist.add_action("strip", "!CA", traj.top)
        alist.add_action("outtraj", "./output/test_trajout.nc", traj.top)
        alist.do_actions([traj[0], traj[1], traj,
                          traj.chunk_iter(chunk=4), traj.frame_iter()])
        Nframes = 2 + traj.n_frames + 8 + traj.n_frames
        dflist.write_all_datafiles()
        print (dslist.size)
        print (dslist[0][:])
        print (dslist[1][:])
        print (dslist.get_dataset(dtype='integer'))
        traj2 = mdio.load("./output/test_trajout.nc", traj.top)
        print (traj.n_frames)
        print ('test_trajout.nc has %s frames' % traj2.n_frames)
        print (traj2[0].n_atoms)
        assert traj2.n_frames == Nframes

if __name__ == "__main__":
    unittest.main()
