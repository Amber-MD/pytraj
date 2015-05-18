from __future__ import print_function
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal, eq
from pytraj.utils.check_and_assert import assert_almost_equal as aa_eq
from pytraj.decorators import no_test, test_if_having
import numpy as np

@test_if_having("MDAnalysis")
def test_load_Amber(self):
    import pytraj.io as io
    from MDAnalysis import Universe

    # load pytraj's TrajectoryIterator
    traj = io.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

    # Fail
    #u = Universe(traj.top.filename, traj.filename, 
    #             format='mdcrd', topology_format='prmtop')
    #print (u)
    #io.load_MDAnalysis(u)

    # OK
    from MDAnalysisTests.datafiles import PSF, DCD
    u = Universe(PSF, DCD)
    print (u)
    io.load_MDAnalysis(u)


if __name__ == "__main__":
    test_load_Amber()
