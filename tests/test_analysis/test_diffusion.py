#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.utils import aa_eq


class TestDiffusion(unittest.TestCase):

    def test_diffusion(self):
        traj = pt.iterload(fn('tz2.ortho.nc'), fn('tz2.ortho.parm7'))

        cm = '''
        parm data/tz2.ortho.parm7
        trajin data/tz2.ortho.nc
        diffusion {0} df {1}
        '''

        options = [
            (':WAT@O nocalc', False),
            (':1000 nocalc', False),
            (':1000 nocalc', True),
            (':1000', False),
            (':1000 nocalc', True),
        ]

        for (mask, individual) in options:
            i_text = 'individual' if individual else ''
            updated_cm = cm.format(mask, i_text)

            state = pt.load_cpptraj_state(updated_cm)
            state.run()

            # make nicer labels
            # label should match to `cm`
            label = 'df'

            for d in state.data[1:]:
                d.key = d.key.replace('[', '').replace(']', '').replace(label,
                                                                        '')

            saved_data = state.data[1:].to_dict()
            data = pt.diffusion(traj, mask, individual=individual).to_dict()

            for key in data.keys():
                if key != 'Label':
                    aa_eq(data[key], saved_data[key])
                else:
                    assert list(data[key]) == list(saved_data[key])


if __name__ == "__main__":
    unittest.main()
