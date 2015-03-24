import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.analysis_dict import AnalysisDict

anadict = AnalysisDict()

class Test(unittest.TestCase):
    def test_0(self):
        # TODO and FIXME: got error:
        # Error: No atoms selected for [*]
        # need to read cpptraj code

        trajin = "./data/md1_prod.Tc5b.x"
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = anadict['rms2d']
        dslist = DataSetList()
        act("trajin ./data/md1_prod.Tc5b.x 2drms rmsout ./output/_test_2drms.dat", traj.top,
            dslist=dslist)
        print (act)

if __name__ == "__main__":
    unittest.main()
