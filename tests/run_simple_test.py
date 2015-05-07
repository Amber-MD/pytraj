import unittest
from pytraj import *
from pytraj.parms import *
from pytraj.trajs import *
from pytraj.datasets import *
from pytraj.common_actions import *

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj import run_tests
        run_tests()

    def test_1(self):
        print ("try to make all action objects")
        from pytraj import adict
        # FIXME, TODO : make failed_list empty
        failed_list = ['createreservoir',]

        for key in adict.keys():
            if key not in failed_list:
                adict[key]

    def test_2(self):
        DataSetList()
        print ("try to make all analysis objects")
        from pytraj import analdict
        failed_list = []

        for key in analdict.keys():
            if key not in failed_list:
                analdict[key]

    def test_3(self):
        print ("try to make all dataset stuff")
        DataSet_double()
        DataSet_float()
        DataSet_integer()
        DataSet_string()
        DataSet_MatrixDbl()
        DataSet_MatrixFlt()
        DataSet_Vector()
        DataSet_Coords()
        DataSet_Coords_REF()
        DataSet_Coords_CRD()
        DataSet_Coords_TRJ()

    def test_4(self):
        print ("try to make structure-related objects")
        Topology()
        Molecule()
        Residue()
        Atom()
        Frame()
        Trajectory()
        TrajectoryIterator()
        TrajinList.TrajinList()

    def test_5(self):
        print ("other stuff. throw all tests don't belong anywhere else here")
        from pytraj import cpptraj_dict
        from pytraj.misc import get_atts
        keys = get_atts(cpptraj_dict)
        cdict = cpptraj_dict.__dict__

        for key in keys:
            if isinstance(cdict[key], dict):
                assert cdict[key].keys() is not None

if __name__ == "__main__":
    unittest.main()
    print ("OK")
    from pytraj.__version__ import __version__
    print (___version__)
