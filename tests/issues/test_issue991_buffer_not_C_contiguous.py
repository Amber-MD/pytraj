#!/usr/bin/env python

from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from pytraj.testing import eq, aa_eq

class Test_Issue_991(unittest.TestCase):
    @unittest.skip('there is no pdb file for Topology yet')
    def test_buffer_not_c_contiguous(self):
        # source code was lightly adapted from jmborr
        # https://github.com/Amber-MD/pytraj/issues/991
        traj = pt.load('data/issue991/short.dcd', 'data/issue991/pdb', 
                       mask='(!:1-256)&(@H1,@H2,@H3,@H4,@H5)')
        
        # Trajectory of the center of mass of the first two residues
        minitraj =  np.empty((2,traj.n_frames,3))
        minitraj[0] = pt.center_of_mass(traj,mask=':1')
        minitraj[1] = pt.center_of_mass(traj,mask=':2')
        minitraj = minitraj.transpose((1,0,2))
        
        # "topology" using the first two atoms
        minitop = traj.top['@1,2']
        
        # Save trajectory
        # make sure there is no ValueError
        pt.write_traj(filename='output/junk.crd', traj=minitraj, top=minitop,
                      overwrite=True)

        # load coord back to make sure we correctly write it
        new_traj = pt.iterload('output/junk.crd', minitop)
        aa_eq(minitraj, new_traj.xyz)


if __name__ == "__main__":
    unittest.main()
