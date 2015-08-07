from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        import numpy as np
        from pytraj.misc import to_amber_mask
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        atm = traj.top("@H=")
        indices = atm.indices
        new_mask = to_amber_mask(indices, mode='int_to_str')
        atm2 = traj.top(new_mask)
        # print(new_mask)
        # print(atm2.indices)
        # print(atm.indices)
        print(atm.indices, atm2.indices)
        assert np.all(atm.indices == atm2.indices)


if __name__ == "__main__":
    unittest.main()
