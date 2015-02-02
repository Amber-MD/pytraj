
import unittest
from time import time
from pytraj.base import *
from pytraj import allactions
from pytraj import io as mdio

class TestActionList(unittest.TestCase):
    def test_run_0(self):
        # load traj
        farray = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
    
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
        alist.add_action(stripact, ArgList("@H*"), toplist, None, dsetlist, dflist)
        alist.add_action(allactions.Action_Distance(), ArgList(":2@CA :3@CA out ./output/test_df.dat"), 
                         toplist, None, dsetlist, dflist)
    
        # 
        print ("test setup_actions")
        print (alist.n_actions)
    
        # do checking
        alist.process(toplist[0])
    
        farray2 = FrameArray()
        frame0 = Frame()
        # testing how fast to do the actions
    
        t0 = time()
        # loop all frames
        # use iterator to make faster loop
        # don't use "for i in range(farray.size)"
        for frame in farray:
            # perform actions for each frame
            frame0 = frame.copy()
            alist.do_actions(frame=frame0)
            #alist.do_actions(frame)
    
            # we need to keep the modified frame in farray2
            #farray2.append(frame)
            farray2.append(frame0)
        print (time() - t0)
    
        # make sure that Action_Strip does its job in stripping
        print (farray2.size)
        assert farray2.n_frames == farray.n_frames
        assert farray2[0].n_atoms != farray[0].n_atoms
    
        # it's time to retrieve the data
        ds = cast_dataset(dsetlist[0], dtype='general')
        print (ds[:10])
        print (dir(dflist))
        dflist.write_all_datafiles()

if __name__ == "__main__":
    unittest.main()
