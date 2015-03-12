import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        # is_ensemble = False
        state = CpptrajState()
        state.toplist.add_parm("./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7")
        state.add_trajin("./data/Test_RemdTraj/rem.nc.000 remdtraj remdtrajtemp 300.0",
                         is_ensemble=False)
        #state.add_trajin("./data/Test_RemdTraj/rem.nc.000")

        state.add_trajout(ArgList("./output/0_test_remd_cpptrajstate.x netcdf"))
        # if is_ensemble == True: save traj for all T
        # else: save traj for only targeted T
        trajlist = state.get_trajinlist()

        # we need to update topology for traj too
        for traj in trajlist:
            traj.top = state.toplist[0]
            print (traj.top)
            print (traj.temperatures)
        state.run()

    def test_1(self):
        # is_ensemble = True
        state = CpptrajState()
        state.toplist.add_parm("./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7")
        state.add_trajin("./data/Test_RemdTraj/rem.nc.000 remdtraj remdtrajtemp 300.0", 
                         is_ensemble=True)
        #state.add_trajin("./data/Test_RemdTraj/rem.nc.000")

        state.add_trajout(ArgList("./output/1_test_remd_cpptrajstate.x netcdf"))
        # if is_ensemble == True: save traj for all T
        # else: save traj for only targeted T
        trajlist = state.get_trajinlist()

        # we need to update topology for traj too
        for traj in trajlist:
            traj.top = state.toplist[0]
            print (traj.top)
            print (traj.temperatures)

        state.write_all_datafiles()
        state.run()

    def test_2(self):
        # is_ensemble = None
        print ("read single traj file")
        state = CpptrajState()
        state.toplist.add_parm("./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7")
        state.add_trajin("./data/Test_RemdTraj/rem.nc.000")
                         
        state.add_trajout(ArgList("./output/2_test_remd_cpptrajstate.x netcdf"))
        # if is_ensemble == True: save traj for all T
        # else: save traj for only targeted T
        trajlist = state.get_trajinlist()

        # we need to update topology for traj too
        for traj in trajlist:
            traj.top = state.toplist[0]
            print (traj.top)
            print (traj.temperatures)

        state.write_all_datafiles()
        state.run()

if __name__ == "__main__":
    unittest.main()
