import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        top = traj.top
        print (list(top.trunc_res_atom_name('@CA')))
        print (list(top.trunc_res_atom_name(0)))
        print (list(top.trunc_res_atom_name('@C')))
        print()
        print (list(top.trunc_res_atom_name(':2-14@C,H,N,O')))

if __name__ == "__main__":
    unittest.main()
