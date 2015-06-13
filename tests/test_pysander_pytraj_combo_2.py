from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import test_if_having 

try:
    import sander
    from parmed.amber.readparm import AmberParm
    has_sander_and_parmed = True
except:
    has_sander_and_parmed = False

class Test(unittest.TestCase):
    @test_if_having("pandas")
    def test_0(self):
        if has_sander_and_parmed:
            from pytraj.misc import get_atts
            import pandas as pd
            from pytraj.dataframe import to_dataframe
            from pytraj.externals.get_pysander_energies import get_pysander_energies

            traj_fn = "./data/md1_prod.Tc5b.x"
            top_fn = "./data/Tc5b.top"
            traj = mdio.iterload(traj_fn, top_fn)

            # parm could be string or AmberParm object
            parm = top_fn

            e_dict = get_pysander_energies(parm=parm, traj=traj, igb=8)
            print (pd.DataFrame(e_dict))
        else:
            print ("require both sander and parmed. Skip test")

if __name__ == "__main__":
    unittest.main()
