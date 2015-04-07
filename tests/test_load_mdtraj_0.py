from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj.utils import Timer
from pytraj import adict
from pytraj import common_actions
from pytraj import io as mdio
from pytraj.utils import has_
from pytraj.misc import get_atts
from pytraj.utils.check_and_assert import assert_almost_equal, eq
from pytraj.six_2 import izip

import numpy as np
from pytraj.load_mdtraj import load_mdtraj

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        if has_("mdtraj") and has_("tables"):
            print ("Testing mdtraj and pytraj")
            import mdtraj as md
            import mdtraj.testing
            traj_filename = mdtraj.testing.get_fn('frame0.h5')
            m_traj = md.load(traj_filename)
            farray = load_mdtraj(m_traj)
            print (farray.top)
            assert isinstance(farray.top, Topology) == True
            print (farray)
            print (farray[0, 0])
            print (farray[0].n_atoms)

            # start assertion
            assert farray.top.n_atoms == m_traj.top.n_atoms
            eq(farray.size, m_traj.n_frames)

            for f_m, f_p in izip(m_traj, farray):
                assert_almost_equal(f_m.xyz.flatten(), f_p.coords)

            with Timer() as t:
                d0 = common_actions.calc_distance("@1 @21", farray)
            print ("time for pytraj_0 = %s" % t.time_gap())

            act = adict['distance']
            dslist = DataSetList()
            with Timer() as t:
                act("@1 @21", farray, dslist=dslist)
            print ("time for pytraj_1= %s" % t.time_gap())

            indices = np.array([[0, 20],])
            with Timer() as t:
                dist_m = md.compute_distances(m_traj, indices, periodic=False)
            print ("time for mdtraj = %s" % t.time_gap())
            print (dist_m)
            assert_almost_equal(d0[:], dist_m[:][0])

            with Timer() as t:
                d0_2 = np.asarray([f.calc_distance(indices) for f in farray]).flatten()
            print ("time for pytraj_2 = %s" % t.time_gap())
            N = 20
            x = d0_2[:N]
            y = d0[:N]
            print (x, y)
            assert_almost_equal(x, y)

        else:
            print ("does not have mdtraj and/or pytables")
            print ("skip test")

if __name__ == "__main__":
    unittest.main()
