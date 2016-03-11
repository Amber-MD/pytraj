from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir


class TestAtomicFluct(unittest.TestCase):

    def test_bfactors(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
        iter_options = {'start': 9, 'stop': 30, 'step': 2}

        bfactors = pt.calc_bfactors(traj(**iter_options))
        s_fname = "/".join((cpptraj_test_dir, "Test_AtomicFluct",
                            "fluct.4.dat.save"))
        saved_bfactors = np.loadtxt(s_fname).T[1]

        aa_eq(saved_bfactors, bfactors.T[1])

        b2 = pt.calc_bfactors(traj(**iter_options), dtype='dataset')
        assert b2[0].key == 'B-factors'

    def test_RMSF(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")

        state = pt.load_batch(traj, '''
        rms first
        average crdset MyAvg
        run
        rms ref MyAvg
        atomicfluct out fluct.agr''')
        state.run()

        t0 = traj[:]
        pt.superpose(t0, ref=0)
        avg = pt.mean_structure(t0)
        pt.superpose(t0, ref=avg)
        data = pt.rmsf(t0)
        aa_eq(data, state.data[-1].values)

    def test_calc_atomicfluct_with_unitcell(self):
        # use iterload for load_batch
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        state = pt.load_cpptraj_state('''
        distance :3 :7
        atomicfluct @CA out output/test.agr
        distance :3 :7''', traj)
        state.run()

        # use `load` method
        t0 = pt.load("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        data = pt.atomicfluct(traj, '@CA')
        aa_eq(data, state.data[-2].values)
        # make sure that traj's coordinates were not altered
        # https://github.com/Amber-MD/pytraj/issues/1166
        aa_eq(pt.distance(t0, ':3 :7'), state.data[-1])
        aa_eq(traj.xyz, t0.xyz)


if __name__ == "__main__":
    unittest.main()
