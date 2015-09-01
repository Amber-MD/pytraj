from __future__ import print_function
import unittest
import numpy as np
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.misc import get_atts


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = DataSetList()
        adict['matrix']("byres @CA", traj, dslist=dslist)
        mat = dslist[0]
        #print(mat)
        #print(get_atts(mat))

        n_residues = traj.top.n_residues
        assert mat.data.shape == (n_residues, n_residues)
        arr0 = np.asarray(mat.data)
        indices = np.where(arr0 == 0.0)[0]
        #print(indices.shape)
        #print(len(indices))

        #print(mat.name)
        #print(mat.legend)


if __name__ == "__main__":
    unittest.main()
