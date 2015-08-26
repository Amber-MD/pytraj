from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj.utils import Timer
from pytraj import adict
from pytraj import common_actions
from pytraj import io
from pytraj.utils import has_
from pytraj.misc import get_atts
from pytraj.utils.check_and_assert import assert_almost_equal as aa_eq
from pytraj.utils.check_and_assert import eq
from pytraj.compat import izip
from pytraj.compat import izip as zip
from pytraj.testing import test_if_having

import numpy as np
from pytraj.io import load_mdtraj


class Test(unittest.TestCase):
    @test_if_having("mdtraj")
    def test_0(self):
        traj = io.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        if has_("mdtraj") and has_("tables"):
            #print("Testing mdtraj and pytraj")
            import mdtraj as md
            import mdtraj.testing
            traj_filename = mdtraj.testing.get_fn('frame0.h5')
            m_traj = md.load(traj_filename)
            farray = load_mdtraj(m_traj)
            #print(farray.top)
            assert isinstance(farray.top, Topology) == True
            # need to set empty Box to match mdtraj's results. ack
            farray.top.box = farray.top.box.__class__()
            #print(farray)
            #print(farray[0, 0])
            #print(farray[0].n_atoms)

            # start assertion
            assert farray.top.n_atoms == m_traj.top.n_atoms
            eq(farray.size, m_traj.n_frames)

            for f_m, f_p in izip(m_traj, farray):
                aa_eq(10 * f_m.xyz.flatten(), f_p.coords)

            with Timer() as t:
                d0 = common_actions.calc_distance(farray, "@1 @21")
            #print("time for pytraj_0 = %s" % t.time_gap())

            act = adict['distance']
            dslist = DataSetList()
            with Timer() as t:
                act("@1 @21", farray, dslist=dslist)
            #print("time for pytraj_1= %s" % t.time_gap())

            indices = np.array([[0, 20], ])
            with Timer() as t:
                dist_m = md.compute_distances(m_traj, indices, periodic=False)
            #print("time for mdtraj = %s" % t.time_gap())
            # convert from "nm" to "angstrom"
            aa_eq(d0[:], 10 * dist_m[:][0])

            with Timer() as t:
                d0_2 = np.asarray([f.calc_distance(indices) for f in farray
                                   ]).flatten()
            #print("time for pytraj_2 = %s" % t.time_gap())
            N = 20
            x = d0_2[:N]
            y = d0[:N]
            ##print (x, y)
            aa_eq(x, y)

        else:
            #print("does not have mdtraj and/or pytables")
            #print("skip test")

    @test_if_having("mdtraj")
    def test_1(self):
        # load water box
        import mdtraj as md
        from mdtraj.testing import get_fn
        m_traj = md.load("./data/tz2.ortho.rst7", top="./data/tz2.ortho.parm7")
        true_traj = io.iterload(
            "./data/tz2.ortho.rst7",
            top="./data/tz2.ortho.parm7")
        traj = io.load_mdtraj(m_traj)
        #print(traj.top.box)
        #print(true_traj.top.box)
        assert traj.n_atoms == m_traj.n_atoms == true_traj.n_atoms

        count = 0
        for f0, f1 in izip(traj, true_traj):
            count += 1
            aa_eq(f0.coords, f1.coords)
            assert f0.box.type == f1.box.type == 'ortho'
        assert count == traj.n_frames == true_traj.n_frames

        traj2 = io.load_mdtraj(m_traj, False)
        #print(traj[0, 0], true_traj[0, 0], traj2[0, 0])
        #print(m_traj.xyz[0, 0])
        aa_eq(traj2.xyz, m_traj.xyz, decimal=3)
        aa_eq(traj.xyz, 10 * m_traj.xyz, decimal=3)

        # provide topology
        traj3 = io.load_mdtraj(m_traj, False, traj2.top)
        aa_eq(traj2.xyz, traj3.xyz)

        # load gro file
        t_gro = md.load(get_fn("frame0.gro"))
        traj = io.load_mdtraj(t_gro, autoconvert=False)
        aa_eq(traj.xyz, t_gro.xyz)


if __name__ == "__main__":
    unittest.main()
