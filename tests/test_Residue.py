import os
import unittest
import pytraj
from pytraj.base import *
from pytraj.core import Residue


class TestResidue(unittest.TestCase):

    def test_0(self):
        top = Topology("./data/Tc5b.top")
        print(top.residuelist[0].first_atom_idx)
        print(top.residuelist[0].last_atom_idx)
        print(top.residuelist[0].n_atoms)

        print(top.residuelist[10].first_atom_idx)
        print(top.residuelist[10].last_atom_idx)
        print(top.residuelist[10].n_atoms)

        print(dir(top.atomlist[0]))
        anames = [atom.__str__() for atom in top.atomlist]
        # print anames

    def test_1(self):
        datadir = "./data/"
        tl = TopologyList()
        tl.add_parm(datadir + "Tc5b.top")
        tl.add_parm(datadir + "HP36.top")
        # tl.List()
        top = tl.get_parm(1)

        print("get 1st residue")
        res1 = next(top.residue_iter())

        # get residue's name
        print(res1)
        # print help(res1)
        # extract residue info
        print("1st atom: %s" % res1.first_atom_idx)
        print("last atom: %s" % res1.last_atom_idx)
        print(res1.original_resnum)
        # print res1.c_str()
        # print res1.Name()
        print(res1.n_atoms)
        # print res1.NameIsSolvent()

        # print all residue's name in top file
        # for res in top.residues():
        #    print res

        # test iterator
        gen = top.residue_iter()
        for idx, res in enumerate(top.residue_iter()):
            print(res.name)
            print(res.n_atoms)
            if idx == 10:
                break
        print(dir(res))
        print(res.first_atom_idx)
        print(res.last_atom_idx)
        print(res.n_atoms)

        res0 = top.residuelist[0]
        print(res0)
        print(top.atomlist[0])
        print(res0.first_atom_idx)
        print(res0.last_atom_idx)
        print(res0.n_atoms)

    def test_4(self):
        top = Topology("./data/Tc5b.top")
        for res in top.residuelist:
            print(res)

    def test_5(self):
        # check segmentation fault for empty Residue
        print(Residue())

if __name__ == '__main__':
    unittest.main()
