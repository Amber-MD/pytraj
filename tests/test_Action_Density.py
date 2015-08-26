from __future__ import print_function
import os
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir

try:
    test_density_dir = os.path.join(cpptraj_test_dir, 'Test_Density')
except:
    test_density_dir = None

delta = 'delta 0.25'
masks = '":PC@P31" ":PC@N31" ":PC@C2" ":PC | :OL | :OL2"'
command = " ".join(["mass out ./output/test_density.dat", delta, masks])


class Test(unittest.TestCase):
    @test_if_path_exists(test_density_dir)
    def test_0(self):
        from pytraj.common_actions import calculate
        # creat mutable Trajectory
        traj = mdio.iterload("./data/DOPC.rst7", "./data/DOPC.parm7")
        farray = traj[:]
        #print(farray)

        # centering
        f0 = farray[0].copy()
        #print(f0.tolist()[:10])
        center = adict['center']
        center('":PC | :OL | :OL2" origin', farray)
        f0 = farray[0].copy()
        #print(f0.tolist()[:10])

        # do action
        act = adict['density']
        dslist = DataSetList()
        act(command, farray, dslist=dslist)
        act.print_output()
        #print(dslist.size)

    @test_if_path_exists(test_density_dir)
    def test_1(self):
        from pytraj import set_world_silent
        set_world_silent(False)
        import pytraj.common_actions as pyca
        from pytraj.common_actions import calculate
        traj = mdio.iterload("./data/DOPC.rst7", "./data/DOPC.parm7")
        fa = traj[:]
        fa.center('":PC | :OL | :OL2" origin')
        command = 'mass delta 0.25 ":PC@P31" ":PC@N31" ":PC@C2" ":PC | :OL | :OL2"'
        dslist = pyca.calc_density(traj, command)
        #print(dslist)


if __name__ == "__main__":
    unittest.main()
