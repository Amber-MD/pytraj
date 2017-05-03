#!/usr/bin/env python
from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from pytraj.testing import aa_eq, tempfolder

# local
from utils import fn


class TestRDF(unittest.TestCase):
    def test_rdf(self):
        traj = pt.iterload(
            fn("tz2.truncoct.nc"),
            fn("tz2.truncoct.parm7"),
            frame_slice=(0, 10))

        command = '''
        radial Radial.agr 0.5 10.0 :5@CD :WAT@O
        radial cRadial.agr 0.5 10.0 :5 :WAT@O center1
        radial cRadial.agr 0.5 10.0 :5 :WAT@O center2
        radial cRadial.agr 0.5 20.0 :3 :WAT@O
        radial cRadial.agr 0.5 20.0 :3 :WAT@O noimage
        radial radial.dat 0.5 10.0 :5@CD :WAT@O
        radial radial2.dat 0.25 10.0 :5@CD :WAT@O
        radial radial2.dat 0.25 10.0 :5@CD :WAT@O volume
        '''

        state = pt.load_batch(traj, command)
        with tempfolder():
            # get data directly from cpptraj
            state.run()

            # get data from pytraj
            data0 = pt.rdf(
                traj,
                solute_mask=':WAT@O',
                bin_spacing=0.5,
                maximum=10.0,
                solvent_mask=':5@CD')

            data01 = pt.rdf(
                traj,
                solvent_mask=':5@CD',
                solute_mask=':WAT@O',
                bin_spacing=0.5,
                maximum=10.0)

            data1 = pt.rdf(
                traj,
                solvent_mask=':5',
                solute_mask=':WAT@O',
                bin_spacing=0.5,
                maximum=10.0,
                center_solvent=True)

            data2 = pt.rdf(
                traj,
                solvent_mask=':5',
                solute_mask=':WAT@O',
                bin_spacing=0.5,
                maximum=10.0,
                center_solute=True)

            data3 = pt.rdf(
                traj,
                solvent_mask=':3',
                solute_mask=':WAT@O',
                bin_spacing=0.5,
                maximum=20.0,
                center_solute=False)

            data4 = pt.rdf(
                traj,
                solvent_mask=':3',
                solute_mask=':WAT@O',
                bin_spacing=0.5,
                maximum=20.0,
                center_solute=False,
                image=False)

            data5 = pt.rdf(
                traj,
                solute_mask=':WAT@O',
                bin_spacing=0.25,
                maximum=10.0,
                solvent_mask=':5@CD')

            # solvent_mask is array
            solvent_indices = pt.select(':WAT@O', traj.top)
            data6 = pt.rdf(
                traj,
                solvent_mask=':5@CD',
                solute_mask=solvent_indices,
                bin_spacing=0.25,
                maximum=10.0)

            # volume
            data7 = pt.rdf(
                traj,
                solvent_mask=':5@CD',
                solute_mask=':WAT@O',
                bin_spacing=0.25,
                maximum=10.0,
                volume=True)

            # do assertion
            aa_eq(data0[1], state.data[1], decimal=7)
            aa_eq(data1[1], state.data[2], decimal=7)
            aa_eq(data2[1], state.data[3], decimal=7)
            aa_eq(data3[1], state.data[4], decimal=7)
            aa_eq(data4[1], state.data[5], decimal=7)
            aa_eq(data7[1], state.data[8], decimal=7)

            # default solvent mask :WAT@O
            aa_eq(data01[1], state.data[1], decimal=7)
            steps = np.loadtxt('radial.dat').T[0]
            aa_eq(data0[0], steps)

            steps2 = np.loadtxt('radial2.dat').T[0]
            aa_eq(data5[0], steps2)
            aa_eq(data6[0], steps2)


if __name__ == "__main__":
    unittest.main()
