import unittest
from time import time
from pytraj.base import *
from pytraj import allactions
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


def run_0():
    # load traj
    farray = mdio.load("../tests/data/md1_prod.Tc5b.x",
                       "../tests/data/Tc5b.top")

    # create 'strip' action
    stripact = allactions.Action_Strip()

    # creat datasetlist to hold distance data
    dsetlist = DataSetList()

    # creat datafilelist to hold filenames and do writing later
    dflist = DataFileList()

    # creat ActionList to hold all actions
    alist = ActionList()

    # creat TopologyList to hold Topology instances
    toplist = TopologyList()
    # add parm
    toplist.add_parm(farray.top)

    # add actions: Action_Strip, Action_Distance and Action_Rmsd
    alist.add_action(stripact, ArgList("@H*"), toplist, dsetlist, dflist)
    alist.add_action(allactions.Action_Distance(),
                     ArgList(":2@CA :3@CA out ./output/_distance.dat"),
                     toplist, dsetlist, dflist)
    alist.add_action(allactions.Action_Rmsd(),
                     ArgList("rms first @CA out ./output/_rmsd.dat"), toplist,
                     dsetlist, dflist)

    #
    print("test setup_actions")
    print("number of actions = ", alist.n_actions)

    # do checking
    alist.process(toplist[0])

    farray2 = Trajectory()

    # creat frame0 as data holder to avoid free memory twices
    # TODO : fix this
    frame0 = Frame()
    # testing how fast to do the actions

    t0 = time()
    # loop all frames
    # use iterator to make faster loop
    # don't use "for i in range(farray.size)"
    for frame in farray:
        # perform actions for each frame
        frame0 = frame.copy()
        alist.do_actions(frame0)
        #alist.do_actions(frame)

        # we need to keep the modified frame in farray2
        #farray2.append(frame)
        farray2.append(frame0)
    print(time() - t0)

    # make sure that Action_Strip does its job in stripping
    print(farray2.size)
    assert farray2.n_frames == farray.n_frames
    assert farray2[0].n_atoms != farray[0].n_atoms

    # it's time to retrieve the data
    # get distance data
    # we need to explicitly cast_dataset
    # future : automatically cast data

    # get distance
    ds0 = cast_dataset(dsetlist[0], dtype='general')
    # get rmsd data
    ds1 = cast_dataset(dsetlist[1], dtype='general')
    print(ds0[:10])
    print(ds1[:10])

    # reproduce cpptraj's output?
    import numpy as np
    rmsdcpp = np.loadtxt("../tests/data/rmsd_to_firstFrame_CA_allres.Tc5b.dat",
                         skiprows=1).transpose()[1][:10]
    # YES
    assert_almost_equal(rmsdcpp, ds1[:10])

    # write output for rmsd and distance (stored in dflist)
    # datatfile: ./_rmsd.dat, _distance.dat
    dflist.write_all_datafiles()
    print(dir(dflist))

    # add more
    # FIXME : "Command terminated" error
    #dflist.add_dataset("./output/dfout_0.dat", ds0)


if __name__ == "__main__":
    run_0()
