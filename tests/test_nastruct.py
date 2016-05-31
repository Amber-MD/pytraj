#!/usr/bin/env python
from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from pytraj.testing import aa_eq, cpptraj_test_dir

class TestNastruct(unittest.TestCase):

    def test_nupars(self):
        fn = "./data/Test_NAstruct/adh026.3.pdb"
        root = 'data/Test_NAstruct/'
        traj = pt.iterload(fn, fn)
        data = pt.nastruct(traj)

        # default
        text = '''
        parm "./data/Test_NAstruct/adh026.3.pdb"
        trajin "./data/Test_NAstruct/adh026.3.pdb"
        nastruct groovecalc 3dna
        '''

        state = pt.load_cpptraj_state(text)
        state.run()

        for key in ['major', 'minor', 'twist']:
            cpp_data = np.array([x.values for x in state.data if x.aspect ==
                                 key])
            # need to transpose to get shape=(n_frames, n_pairs)
            cpp_data = cpp_data.T
            aa_eq(data[key][1], cpp_data)

        # TODO: assert
        data._summary(np.mean, indices=None)
        data._summary(np.mean, indices=[1, ])
        data._summary(np.mean, keys=['major', 'twist'], indices=[1, ])
        data._summary(np.mean, keys='major', indices=[1, ])
        data._summary(np.std, indices=[1, ])
        data._summary([np.std, np.mean], indices=[1, ])
        data._explain()
        dir(data)

        # pickle
        pt.to_pickle(data, 'data/na.pk')
        na2 = pt.read_pickle('data/na.pk')

        for key in data.keys():
            aa_eq(data[key][1], na2[key][1])

        # raise
        self.assertRaises(ValueError, lambda: pt.nastruct(traj, dtype='ndarray'))

    def test_baseref(self):
        fn = "./data/Test_NAstruct/adh026.3.pdb"
        baseref_fn = 'data/Test_NAstruct/Atomic_G.pdb.nastruct'

        traj = pt.iterload(fn, fn)
        data = pt.nastruct(traj, baseref=baseref_fn)

        saved_data = """
         # major minor 3 frames
         19.9015, 15.8085,
         20.7510, 15.8180,
         19.9822, 15.4856,

         19.5528, 15.9407,
         20.6965, 16.2324,
         19.6335, 16.3369,

         19.2871, 15.9514,
         21.2010, 16.7085,
         19.9055, 16.6283
        """

        aa_eq(data.major[1][0], [19.9015, 20.7510, 19.9822], decimal=3)
        aa_eq(data.major[1][1], [19.5528, 20.6965, 19.6335], decimal=3)
        aa_eq(data.major[1][2], [19.2871, 21.2010, 19.9055], decimal=3)

        aa_eq(data.minor[1][0], [15.8085, 15.8180, 15.4856], decimal=3)

    def test_nupars_vs_x3dna(self):
        traj = pt.iterload('data/Test_NAstruct/x3dna/rna.pdb')
        ref = pt.iterload('data/Test_NAstruct/x3dna/rna_nab.pdb')
        nu = pt.nastruct(traj, ref=ref, groove_3dna=True)

        root = 'data/Test_NAstruct/x3dna/'

        # helical pars
        saved_helical_pars = np.loadtxt(root + 'bp_helical.par',
                                        skiprows=3,
                                        usecols=range(1, 13)).T
        aa_eq(nu.shear[1], saved_helical_pars[0], decimal=3)
        aa_eq(nu.stretch[1], saved_helical_pars[1], decimal=3)
        aa_eq(nu.stagger[1], saved_helical_pars[2], decimal=3)
        aa_eq(nu.buckle[1], saved_helical_pars[3], decimal=3)
        aa_eq(nu.prop[1], saved_helical_pars[4], decimal=3)
        aa_eq(nu.open[1], saved_helical_pars[5], decimal=3)
        aa_eq(nu.xdisp[1], saved_helical_pars[6][1:], decimal=3)
        aa_eq(nu.ydisp[1], saved_helical_pars[7][1:], decimal=3)
        aa_eq(nu.hrise[1], saved_helical_pars[8][1:], decimal=3)
        aa_eq(nu.incl[1], saved_helical_pars[9][1:], decimal=3)
        aa_eq(nu.tip[1], saved_helical_pars[10][1:], decimal=3)
        aa_eq(nu.htwist[1], saved_helical_pars[11][1:], decimal=3)

        # bp_step
        saved_helical_pars = np.loadtxt(root + 'bp_step.par',
                                        skiprows=3,
                                        usecols=range(1, 13)).T
        aa_eq(nu.shift[1], saved_helical_pars[6][1:], decimal=3)
        aa_eq(nu.slide[1], saved_helical_pars[7][1:], decimal=3)
        aa_eq(nu.tilt[1], saved_helical_pars[9][1:], decimal=3)
        aa_eq(nu.roll[1], saved_helical_pars[10][1:], decimal=3)
        aa_eq(nu.twist[1], saved_helical_pars[11][1:], decimal=3)

        # grove
        aa_eq(nu.minor[1][0], [15.8, 15.8, 15.5], decimal=1)
        aa_eq(nu.major[1][0], [19.9, 20.8, 20.0], decimal=1)


if __name__ == "__main__":
    unittest.main()
