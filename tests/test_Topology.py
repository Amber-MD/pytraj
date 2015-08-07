import os
import numpy as np
import unittest
from pytraj.Topology import Topology
from pytraj.AtomMask import AtomMask
from pytraj.core.FileName import FileName
from pytraj.base import *
from pytraj.decorators import no_test
from pytraj.core.NameType import NameType

TRAJ = Trajectory("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")


class TestTopology(unittest.TestCase):
    def test_empty_top(self):
        top = Topology()
        print(top.file_path())
        assert top.is_empty() == True
        filename = "./data/Tc5b.top"
        top = Topology(filename)
        print(top.file_path())
        assert top.is_empty() == False

    def test_1(self):
        datadir = "./data/"
        filename = "./data/Tc5b.top"

        top = Topology(filename)
        #top2 = top.modify_state_by_mask(AtomMask("!@CA"))
        #
        top.summary()
        top.atom_info("@CA")

        print("test strip_atoms: strip all but CA")
        top.strip_atoms("!@CA")
        assert top.n_atoms == 20
        top.summary()
        top.atom_info("@CA")
        top.atom_info("@H")

        print("test Parm_Amber for write")
        #Parm_Amber().write_parm("test_write.top", top)
        # top.write_parm("test_write2.top")

        print("test atom_iterator")
        for atom in top.atom_iter():
            print(atom)

        print("test res_iterator")
        for res in top.residue_iter():
            pass

        print("test mol_iterator")
        for mol in top.mol_iter():
            pass

        for idx, atom in enumerate(top.atoms):
            print(atom)
        assert idx + 1 == top.n_atoms

    def test_select_mask(self):
        top = Topology("./data/Tc5b.top")
        arr0 = top.atom_indices("@CA")
        print(arr0)
        print(type(arr0))

    def test_call(self):
        traj = Trajectory("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        #top = traj.top
        #frame = traj[0]
        # print(frame[top(":2-18@CA")])

    def test_get_unique(self):
        traj = Trajectory("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        top = traj.top
        print(top.atom_names)
        print(top.residue_names)

    def test_len(self):
        traj = Trajectory("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        top = traj.top
        assert len(top) == top.n_atoms

    def test_charge(self):
        print(TRAJ.top.total_charge)


if __name__ == "__main__":
    unittest.main()
