import os
import unittest
from pytraj.Topology import Topology
from pytraj.ArgList import ArgList
from pytraj.Frame import Frame
from pytraj.AtomMask import AtomMask
from pytraj.FrameList import FrameList
from pytraj.TopologyList import TopologyList 

class TestFrameList(unittest.TestCase):
    def test_1(self):
        datadir = "./data/"
        toplist = TopologyList()
        toplist.add_parm("./data/Tc5b.top")
        top = toplist.get_parm(0)
        arglist = ArgList("./data/Tc5b.nat.crd")
        print("====================")
        print("toplist.info()")
        toplist.info()
        print("====================")
        try:
            top.n_atoms = 3000
        except:
            print("can't not set n_atoms")
        print("====================")
        
        #mylist = FrameList()
        #mylist.add_ref_frame(arglist, toplist)
        #aref = mylist.get_active_reference()
        #mylist.list()
        #print mylist.n_frames
        #
        #print aref.xyz

if __name__ == "__main__":
    unittest.main()
