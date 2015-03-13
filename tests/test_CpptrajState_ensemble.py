import unittest
from itertools import izip
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    # TODO : need to change `is_ensemble` to something else
    # what does it mean with "remdtraj" and "is_ensemble = False" together
    def test_0(self):
        # is_ensemble = False
        print ("read esemble for single temperature")
        state = CpptrajState()
        state.toplist.add_parm("./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7")
        state.add_trajin("./data/Test_RemdTraj/rem.nc.000 remdtraj remdtrajtemp 300.0",
                         is_ensemble=False)
        #state.add_trajin("./data/Test_RemdTraj/rem.nc.000")

        state.add_trajout(ArgList("./output/0_test_remd_cpptrajstate.x netcdf"))
        # if is_ensemble == True: save traj for all T
        # else: save traj for only targeted T
        trajlist = state.get_trajinlist()
        print (trajlist.size)

        # we need to update topology for traj too
        for traj in trajlist:
            traj.top = state.toplist[0]
            print (traj.top)
            print (traj.temperatures)
        state.run()

    def test_1(self):
        # is_ensemble = True
        print ("read esemble for all temperatures")
        state = CpptrajState()
        state.toplist.add_parm("./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7")
        state.add_trajin("./data/Test_RemdTraj/rem.nc.000 remdtraj remdtrajtemp 300.0", 
                         is_ensemble=True)
        #state.add_trajin("./data/Test_RemdTraj/rem.nc.000")

        state.add_trajout(ArgList("./output/1_test_remd_cpptrajstate.x netcdf"))
        # if is_ensemble == True: save traj for all T
        # else: save traj for only targeted T
        trajlist = state.get_trajinlist()
        print (trajlist.size)

        traj0 = trajlist[0]
        traj0.top = state.toplist[0]

        saved_traj = mdio.load("data/Test_RemdTraj/temp0.crd.300.00", 
                               "./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7")

        # make sure that we DO get 300K traj
        for f0, f1 in izip(traj0, saved_traj):
            print (f0, f1)
            assert_almost_equal(f0.coords, f1.coords)

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
        print (trajlist.size)

        # we need to update topology for traj too
        for traj in trajlist:
            traj.top = state.toplist[0]
            print (traj.top)
            print (traj.temperatures)

        state.write_all_datafiles()
        state.run()

if __name__ == "__main__":
    unittest.main()
