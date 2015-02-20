import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.Energy import Energy_Amber

class Test(unittest.TestCase):
    def test_0(self):
        # FIXME: correct AtomMask object. currently get '0' for energies
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        ene = Energy_Amber()
        atm = AtomMask("@CA")
        traj.top.set_char_mask(atm)
        print (atm.n_atoms)
        print (atm.selected_indices())
        print (ene.E_bond(traj[0], traj.top, atm))
        print (ene.E_angle(traj[0], traj.top, atm))

if __name__ == "__main__":
    unittest.main()
