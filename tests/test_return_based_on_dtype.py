from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
from pytraj.utils.check_and_assert import is_word_in_class_name as is_in_class_name
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    def test_0(self):
        import numpy as np
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        d0 = traj.calc_COG(dtype='dataset') # default = 'dataset'
        print (d0)
        assert is_in_class_name(d0, 'DatasetVector') == True
        mylist = traj.calc_COG(dtype='list') # default = 'dataset'
        print (mylist)
        assert isinstance(mylist, list) == True
        mynp = traj.calc_COG(dtype='ndarray') # default = 'dataset'
        assert is_in_class_name(mynp, 'ndarray')
        aa_eq(np.asarray(mylist).flatten(), mynp.flatten())

if __name__ == "__main__":
    unittest.main()
