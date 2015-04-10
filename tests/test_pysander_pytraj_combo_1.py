from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import test_if_having 

try:
    import sander
    from chemistry.amber.readparm import AmberParm
    has_sander_and_parmed = True
except:
    has_sander_and_parmed = False

class Test(unittest.TestCase):
    @test_if_having("pandas")
    def test_0(self):
        if has_sander_and_parmed:
            from pytraj.misc import get_atts
            from collections import defaultdict
            ddict = defaultdict(list, [])

            traj_fn = "./data/md1_prod.Tc5b.x"
            top_fn = "./data/Tc5b.top"
            traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
            parm = AmberParm(top_fn)
            inp = sander.gas_input(8)
            parm.load_coordinates(traj[0].coords)

            with sander.setup(parm, parm.coords, None, inp):
                for frame in traj:
                    sander.set_positions(frame.coords)
                    #sander.set_positions(frame.buffer1d)
                    ene, frc = sander.energy_forces()
                    ene_atts = get_atts(ene)
                    for att in ene_atts:
                        ddict[att].append(getattr(ene, att))
            assert sander.is_setup() == False

            # cpptraj
            from pytraj import DataSetList
            dslist = DataSetList()
            act = adict['energy']
            act("", traj, dslist=dslist)
            print (dslist['ENE_00000[bond]'][0][:])

            # to DataFrame
            import pandas as pd
            from pytraj.dataframe import to_dataframe
            dframe = to_dataframe(dslist)
            print (dframe)
            print (pd.DataFrame(ddict))

        else:
            print ("require both sander and parmed. Skip test")


if __name__ == "__main__":
    unittest.main()
