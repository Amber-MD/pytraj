#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from numpy.testing import assert_equal


class TestDiffusion(unittest.TestCase):

    def test_diffusion(self):
        traj = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')

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
