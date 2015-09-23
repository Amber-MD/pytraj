#!/usr/bin/env python
from __future__ import print_function
import unittest
import numpy as np
import os
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir


class TestRDF(unittest.TestCase):
    def test_rdf(self):
        traj = pt.iterload(
            "./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7", frame_slice=(0, 10))

        command = '''
        radial Radial.agr 0.5 10.0 :5@CD :WAT@O
        radial cRadial.agr 0.5 10.0 :5 :WAT@O center1
        radial cRadial.agr 0.5 10.0 :5 :WAT@O center2
        radial cRadial.agr 0.5 20.0 :3 :WAT@O
        radial cRadial.agr 0.5 20.0 :3 :WAT@O noimage
        '''

        # get data directly from cpptraj
        state = pt.load_batch(traj, command)
        state.run()
        
        # get data from pytraj
        data0 = pt.rdf(traj, solvent_mask=':WAT@O', bin_spacing=0.5,
                      maximum=10.0,
                      solute_mask=':5@CD')

        data01 = pt.rdf(traj, solvent_mask=':WAT@O', bin_spacing=0.5,
                      maximum=10.0,
                      solute_mask=':5@CD')

        data1 = pt.rdf(traj, solvent_mask=':WAT@O', bin_spacing=0.5,
                      maximum=10.0,
                      center_solvent=True,
                      solute_mask=':5')

        data2 = pt.rdf(traj, solvent_mask=':WAT@O', bin_spacing=0.5,
                      maximum=10.0,
                      center_solute=True,
                      solute_mask=':5')

        data3 = pt.rdf(traj, solvent_mask=':WAT@O', bin_spacing=0.5,
                      maximum=20.0,
                      center_solute=False,
                      solute_mask=':3')

        data4 = pt.rdf(traj, solvent_mask=':WAT@O', bin_spacing=0.5,
                      maximum=20.0,
                      center_solute=False,
                      image=False,
                      solute_mask=':3')

        # do assertion
        aa_eq(data0, state.data[0], decimal=7)
        aa_eq(data1, state.data[1], decimal=7)
        aa_eq(data2, state.data[2], decimal=7)
        aa_eq(data3, state.data[3], decimal=7)
        aa_eq(data4, state.data[4], decimal=7)

        # default solvent mask :WAT@O
        aa_eq(data01, state.data[0], decimal=7)

if __name__ == "__main__":
    unittest.main()
