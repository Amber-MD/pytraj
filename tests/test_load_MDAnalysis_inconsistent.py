from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal, eq
from pytraj.utils.check_and_assert import assert_almost_equal as aa_eq
from pytraj.decorators import no_test, test_if_having
import numpy as np


@test_if_having("MDAnalysis")
def test_load_Amber():
    import numpy as np
    import pytraj.io as io
    from MDAnalysis import Universe

    # load pytraj's TrajectoryIterator
    traj = io.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

    u = Universe(traj.top.filename, traj.filename,
                 format='mdcrd',
                 topology_format='prmtop')
    utraj = io.load_MDAnalysis(u)
    aa_eq(np.asarray(utraj.calc_dssp()).flatten(),
          np.asarray(traj.calc_dssp()).flatten())

    # OK
    from MDAnalysisTests.datafiles import PSF, DCD
    u = Universe(PSF, DCD)
    print(u)
    io.load_MDAnalysis(u)


if __name__ == "__main__":
    test_load_Amber()
