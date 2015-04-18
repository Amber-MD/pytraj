import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
top = traj.top

class Test(unittest.TestCase):
    def test_0(self):
        print (list(top.trunc_res_atom_name('@CA')))
        print (list(top.trunc_res_atom_name(0)))
        print (list(top.trunc_res_atom_name('@C')))
        print()
        print (list(top.trunc_res_atom_name(':2-14@C,H,N,O')))
        
    def test_1(self):
        print ("find_atom_in_residue")
        name = "CA  "
        print (top.find_atom_in_residue(3, name))
        print (top[58].name == name)

    def test_2(self):
        print ("test top[indices]")
        indices = top("@CA").indices
        indices_s = top.select("@CA")
        assert indices_s == indices
        atom_list = top[indices]

        assert len(atom_list) == len(indices)

        for atom in atom_list:
            assert atom.name == 'CA  '

if __name__ == "__main__":
    unittest.main()
