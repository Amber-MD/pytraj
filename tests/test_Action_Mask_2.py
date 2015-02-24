import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj.actions.Action_Mask_2 import Action_Mask_2
        traj = mdio.load("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        act = Action_Mask_2()
        dslist = DataSetList()
        dflist = DataFileList()

        maskout = "./output/tz2_Action_Mask.mask.out"
        maskpdb = "./output/tz2_Action_Mask.mask.pdb"
        command = 'mask "(:5 <:3.0) & :WAT" maskout %s maskpdb %s' % (maskout, maskpdb)
        command2 = 'mask "(:5 <:3.0) & :WAT" maskout mask.out maskmol2 mask.mol2'

        act(command,
            traj(0, 0), traj.top, dslist=dslist,
            dflist=dflist)
        print (dslist.size)
        print (dslist[0][:])
        print (dslist[1][:])
        print (dslist[0].name)
        print (dslist[1].name)

        #dflist.write_all_datafiles()

        #act(command2,
        #    traj(), traj.top,
        #    dflist=dflist)

        #dflist.write_all_datafiles()

if __name__ == "__main__":
    unittest.main()
