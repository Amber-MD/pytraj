import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.analyses.Analysis_Rms2d import Analysis_Rms2d

class Test(unittest.TestCase):
    def test_0(self):
        top = Topology("./data/Tc5b.top")
        arg = """
        parm ./data/Tc5b.top
        trajin ./data/md1_prod.Tc5b.x
        @CA rmsout ./output/_test_rms2d.dat
        """
        act = Analysis_Rms2d()
        act(arg, top)

if __name__ == "__main__":
    unittest.main()
