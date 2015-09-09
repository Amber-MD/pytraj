from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

try:
    import sander
    from parmed.amber.readparm import AmberParm
    has_sander_and_parmed = True
except:
    has_sander_and_parmed = False


class Test(unittest.TestCase):
    def test_0(self):
        if has_sander_and_parmed:
            traj_fn = "./data/md1_prod.Tc5b.x"
            top_fn = "./data/Tc5b.top"
            traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
            parm = AmberParm(top_fn)
            inp = sander.gas_input(8)

            for frame in traj:
                parm.load_coordinates(frame.coords)
                sander.setup(parm, parm.coordinates, None, inp)
                ene, frc = sander.energy_forces()
                sander.cleanup()
        else:
            pass

    def test_1(self):
        if has_sander_and_parmed:
            traj_fn = "./data/md1_prod.Tc5b.x"
            top_fn = "./data/Tc5b.top"
            traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
            parm = AmberParm(top_fn)
            inp = sander.gas_input(8)

            import numpy as np

            for frame in traj:
                arr0 = np.asarray(frame.buffer1d)
                parm.load_coordinates(arr0)
                sander.setup(parm, parm.coordinates, None, inp)
                ene, frc = sander.energy_forces()
                sander.cleanup()
        else:
            pass


if __name__ == "__main__":
    unittest.main()
