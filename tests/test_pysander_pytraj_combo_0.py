from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

import sander
from chemistry.amber.readparm import AmberParm

class Test(unittest.TestCase):
    def test_0(self):
        traj_fn = "./data/md1_prod.Tc5b.x"
        top_fn = "./data/Tc5b.top"
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        parm = AmberParm(top_fn)
        inp = sander.gas_input(8)

        for frame in traj:
            parm.load_coordinates(frame.coords)
            sander.setup(parm, parm.coords, None, inp)
            ene, frc = sander.energy_forces()
            print (ene.gb)
            sander.cleanup()

if __name__ == "__main__":
    unittest.main()
