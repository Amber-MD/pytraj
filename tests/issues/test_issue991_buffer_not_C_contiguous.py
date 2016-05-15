#!/usr/bin/env python

from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from pytraj.testing import eq, aa_eq


class TestIssue991(unittest.TestCase):

    def test_buffer_not_c_contiguous(self):
        # source code was lightly adapted from jmborr
        # https://github.com/Amber-MD/pytraj/issues/991
        traj = pt.load('data/issue991/short.dcd',
                       'data/issue991/pdb.gz',
                       mask='(!:1-256)&(@H1,@H2,@H3,@H4,@H5)')

        # Trajectory of the center of mass of the first two residues
        minitraj = np.empty((2, traj.n_frames, 3))
        minitraj[0] = pt.center_of_mass(traj, mask=':1')
        minitraj[1] = pt.center_of_mass(traj, mask=':2')
        minitraj = minitraj.transpose((1, 0, 2))

        # "topology" using the first two atoms
        minitop = traj.top['@1,2']

        # Save trajectory
        # make sure there is no ValueError
        # something is wrong with pdb, crd extension when loading with
        # minitop (poor topology?)
        # make issue in cpptraj too?

        exts = ['nc', 'dcd', 'mdcrd', 'crd', 'pdb', 'trr']
        exts.remove('mdcrd')
        exts.remove('crd')

        for ext in exts:
            fn = 'output/junk.' + ext
            pt.write_traj(filename=fn,
                          traj=minitraj,
                          top=minitop,
                          overwrite=True)

            # load coord back to make sure we correctly write it
            new_traj = pt.iterload(fn, minitop)
            # mdcrd, crd, pdb has only 3 digits after decimal
            decimal = 5 if ext in ['nc', 'dcd', 'trr'] else 3
            aa_eq(minitraj, new_traj.xyz, decimal=decimal)


if __name__ == "__main__":
    unittest.main()
