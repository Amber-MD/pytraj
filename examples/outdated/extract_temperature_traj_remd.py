"""This script shows how to extract frames having the same temperature
from replica exchange MD run. You can do it with cpptraj but this shows how easily
to write new script with pytraj
# TODO : check typos for DOCo
"""

import unittest
from array import array
from glob import glob
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


def get_frames_same_T():
    # create a list of all remd trajs
    flist = glob("../tests/data/Test_RemdTraj/rem.nc.*")

    # make a list of TrajReadOnly instances
    trajlist = []
    for fh in flist:
        topfile = "../tests/data/Test_RemdTraj/ala2.99sb.mbondi2.parm7"

        # load trajectory and append to trajlist
        trajlist.append(mdio.load(fh, topfile))

    # make Trajectory instance that holds 492.2 T frames
    # we need to reserve n_frames to hold the data
    f4922 = Trajectory(n_frames=trajlist[0].n_frames)

    assert f4922.n_frames == trajlist[0].n_frames
    f4922.top = trajlist[0].top.copy()

    # extract frames having T = 492.2
    # use iteration for nested loops
    for traj in trajlist:
        for idx, frame in enumerate(traj):
            if frame.temperature == 492.2:
                # we don't use `append` method since we want to make sure
                # frames are in the order of simulation time
                f4922[idx] = frame

    # make sure f4922 only hold frames having T = 492.2 K
    arr0 = array('d', [492.2, 492.2, 492.2, 492.2, 492.2, 492.2, 492.2, 492.2,
                       492.2, 492.2])
    assert f4922.temperatures == arr0

    # make sure we reproduce cpptraj output
    cpptraj = mdio.load("../tests/data/Test_RemdTraj/temp0.crd.492.20",
                        topfile)
    print(f4922[5].coords[:10])
    print(cpptraj[5].coords[:10])
    for idx, framepy in enumerate(f4922):
        assert_almost_equal(framepy.coords, cpptraj[idx].coords)
        print("rmsd between pytraj Frame and cpptraj Frame = %s " %
              framepy.rmsd(cpptraj[idx]))

    # FIXME: `rmsd` do the fitting in the fly
    # coords of frame would be changed
    print("YES, we can reproduce cpptraj output")


if __name__ == "__main__":
    get_frames_same_T()
