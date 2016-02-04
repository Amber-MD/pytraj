from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.utils.check_and_assert import is_word_in_class_name
from pytraj.utils import has_


class Test(unittest.TestCase):

    def test_0(self):
        from pytraj.all_actions import calc_matrix
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        d0 = calc_matrix(traj, "@CA", dtype='dataset')
        assert is_word_in_class_name(d0, 'DatasetList') == True


if __name__ == "__main__":
    unittest.main()
