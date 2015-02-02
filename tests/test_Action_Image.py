import unittest
from pytraj.base import *
from pytraj import allactions
from pytraj.cast_dataset import cast_dataset
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.utils.check_and_assert import file_exist


# TODO : not reproduce cpptraj's output yet
class TestRadgyr(unittest.TestCase):
    def test_0(self):
        farray = FrameArray(filename="./data/tz2.truncoct.nc", top="./data/tz2.truncoct.parm7")[:2]
        fold = farray.copy()
        print("old file: ", fold[0, 0, :])

        act = allactions.Action_Image()
        ptrajin = """
        center :2-11
        image center familiar com :6
        """

        act.read_input(ptrajin, farray.top)
        act.process(farray.top)
       
        f2 = FrameArray()
        f2.top = farray.top.copy()

        act.do_action(0, farray)

        print(farray[0, 0, :])
        print(fold[0, 0, :])

        if file_exist("./CpptrajTest/Test_Image/image4.crd"):
            fnew = mdio.load("./CpptrajTest/Test_Image/image4.crd", 
                             "./data/tz2.truncoct.parm7")
            assert fnew.size == 2
            # make sure that pytraj reproduce cpptraj outputo
            # TODO : not yet. double-check
            print(fnew[0].same_coords_as(farray[0]))
            print(fnew[0, 0, :])
            print(fold[0].same_coords_as(farray[0]))
            print(fnew[0].rmsd(farray[0]))

if __name__ == "__main__":
    unittest.main()
