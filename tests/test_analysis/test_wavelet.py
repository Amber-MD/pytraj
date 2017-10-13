#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq


class TestWavelet(unittest.TestCase):
    def test_wavelet(self):
        traj = pt.load(fn('DPDP.nc'), fn('DPDP.parm7'))
        traj.superpose('@C,CA,N', ref=0)
        command = 'name mywave nb 10 s0 2 ds 0.25 type morlet correction 1.01 chival 0.25 :1-22 cluster minpoints 66 epsilon 10.0'
        data = pt.wavelet(traj, command)

        state = pt.load_cpptraj_state("""
        parm {}
        trajin {}
        rms @C,CA,N first
        wavelet name mywave nb 10 s0 2 ds 0.25 type morlet correction 1.01 chival 0.25 :1-22 cluster minpoints 66 epsilon 10.0
        """.format(fn('DPDP.parm7'), fn('DPDP.nc')))

        state.run()

        # ignore Topology, COORDS dataset
        cpp_data = state.data[3:].to_dict()
        assert sorted(cpp_data.keys()) == sorted(data.keys())

        for (k0, v0) in data.items():
            aa_eq(v0, cpp_data[k0])


if __name__ == "__main__":
    unittest.main()
