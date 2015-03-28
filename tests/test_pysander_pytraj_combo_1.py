from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

try:
    import sander
    from chemistry.amber.readparm import AmberParm
    has_sander_and_parmed = True
except:
    has_sander_and_parmed = False

class Test(unittest.TestCase):
    def test_0(self):
        if has_sander_and_parmed:
            traj_fn = "./data/md1_prod.Tc5b.x"
            top_fn = "./data/Tc5b.top"
            traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
            parm = AmberParm(top_fn)
            inp = sander.gas_input(8)
            parm.load_coordinates(traj[0].coords)

            with sander.setup(parm, parm.coords, None, inp):
                for frame in traj:
                    #sander.set_positions(frame.coords)
                    sander.set_positions(frame.buffer1d)
                    ene, frc = sander.energy_forces()
                    print (ene.gb)
            assert sander.is_setup() == False
        else:
            print ("require both sander and parmed. Skip test")


if __name__ == "__main__":
    unittest.main()
