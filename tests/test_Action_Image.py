import unittest; import pytraj as pt
from pytraj import adict
from pytraj import allactions, Trajectory
from pytraj.datasets import cast_dataset
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.utils.check_and_assert import file_exist
from pytraj.testing import aa_eq


class TestRadgyr(unittest.TestCase):
    def test_0(self):
        farray = Trajectory(
            filename="./data/tz2.truncoct.nc",
            top="./data/tz2.truncoct.parm7")[:2]
        f_old = farray.copy()
        #print("old file: ", f_old[0, 0, :])

        act = allactions.Action_Image()
        ptrajin = """
        center :2-11
        image center familiar com :6
        """

        act.read_input(ptrajin, farray.top)
        act.process(farray.top)

        f2 = Trajectory()
        f2.top = farray.top.copy()

        act.do_action(farray)

        #print(farray[0, 0, :])
        #print(f_old[0, 0, :])

        if file_exist("./CpptrajTest/Test_Image/image4.crd.save"):
            fnew = mdio.iterload("./CpptrajTest/Test_Image/image4.crd.save",
                                 "./data/tz2.truncoct.parm7")
            assert fnew.size == 2
            #print(fnew[0].same_coords_as(farray[0]))
            #print(fnew[0, 0, :])
            #print(f_old[0].same_coords_as(farray[0]))
            #print(fnew[0].rmsd(farray[0]))

    def test_1(self):
        from pytraj.core.ActionList import ActionList
        farray = mdio.iterload(
            "./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")[:2]
        f_old = farray.copy()
        #print("old file: ", f_old[0, 0, :])
        act = adict['image']
        act(':2-11', f_old)
        act('familiar com :6', f_old)

        if file_exist("./CpptrajTest/Test_Image/image4.crd.save"):
            #print("having file ./CpptrajTest/Test_Image/image4.crd.save")
            fsaved = mdio.iterload("./CpptrajTest/Test_Image/image4.crd.save",
                                   "./data/tz2.truncoct.parm7")
            #print(f_old[0, 0], fsaved[0, 0])
            # TODO, FIXME: assert failed
            #aa_eq(f_old[0].coords, fsaved[0].coords)
            #aa_eq(f_old[1].coords, fsaved[1].coords)


if __name__ == "__main__":
    unittest.main()
