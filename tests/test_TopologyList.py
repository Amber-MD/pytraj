import os
import unittest
from pytraj.TopologyList import TopologyList 
from pytraj.Topology import Topology

class TestTopologyList(unittest.TestCase):
    def test_1(self):
        datadir = "./data/"
        tlist = TopologyList()

        # create TopologyList instance with 4 Topology instances
        tlist.add_parm(Topology())
        tlist.add_parm(datadir + "Tc5b.top")
        tlist.add_parm(datadir + "HP36.top")
        tlist.add_parm(Topology())
        tlist.info()
        
        # test assignment
        tlist[0] = tlist[1]
        assert tlist[0].n_atoms == tlist[1].n_atoms
        print(tlist[0].n_atoms)
        #tlist.info()
        
        #tlist[3].summary()
        
        t_0 = tlist[0]
        t_0.copy(tlist[2])
        #t_0.summary()

        # make sure changing t_0 does not affect tlist[0] since we make a copy
        assert t_0.n_atoms != tlist[0].n_atoms
        #t_0.summary()
        
        t_0_copy = t_0.copy()
        assert t_0_copy != t_0
        assert t_0_copy.n_atoms == t_0.n_atoms

        print("make sure changing t_0_copy does not affect t_0")
        t_0_copy = Topology()
        assert t_0_copy.n_atoms != t_0.n_atoms
        #t_0.summary()
        #t_0_copy.summary()

        # clear tlist
        tlist.clear()
        # testing index out of range
        #tlist[0]
        tlist.info()

        print("adding new Topology instances")
        tlist.add_parm(Topology())
        tlist.add_parm(Topology())
        tlist.add_parm(Topology())
        tlist[0]
        tlist[1]
        tlist[2]

        print("test iteration")
        print("test adding top instances from python list")
        tmplist = [Topology(), Topology(), Topology()]
        tlist_2 = TopologyList()
        tlist_2.get_parm_from_pylist(tmplist)

if __name__ == "__main__":
    unittest.main()
