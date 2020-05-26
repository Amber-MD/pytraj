from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
import pytest


class TestWatershell(unittest.TestCase):
    def test_watershell(self):
        traj = pt.iterload(fn('tz2.truncoct.nc'), fn('tz2.truncoct.parm7'))
        state = pt.load_batch(traj, '''
        watershell :1-7
        ''')

        d0 = pt.watershell(traj, solute_mask=':1-7')
        state.run()
        aa_eq(d0.values, state.data[[1, 2]].values)

        # need to provide solute_mask
        with pytest.raises(ValueError):
            pt.watershell(traj)


if __name__ == "__main__":
    unittest.main()
