import os
import unittest
from pytraj.core.TopologyList import TopologyList
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
        # test assignment
        tlist[0] = tlist[1]
        assert tlist[0].n_atoms == tlist[1].n_atoms
        # tlist.info()

        # tlist[3].summary()

        t_0 = tlist[0]
        t_0 = tlist[2].copy()
        # t_0.summary()

        # make sure changing t_0 does not affect tlist[0] since we make a copy
        assert t_0.n_atoms != tlist[0].n_atoms
        # t_0.summary()

        t_0_copy = t_0.copy()
        assert t_0_copy != t_0
        assert t_0_copy.n_atoms == t_0.n_atoms

        t_0_copy = Topology()
        assert t_0_copy.n_atoms != t_0.n_atoms
        # t_0.summary()
        # t_0_copy.summary()

        tlist.add_parm(Topology())
        tlist.add_parm(Topology())
        tlist.add_parm(Topology())
        tlist[0]
        tlist[1]
        tlist[2]

        tmplist = [Topology(), Topology(), Topology()]
        tlist_2 = TopologyList()


if __name__ == "__main__":
    unittest.main()
