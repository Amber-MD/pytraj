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
        # Aim: `cpptraj` use `first` frame as default for some Actions (Action_Rmsd,
        # Action_NativeContacts, ...). We can specify reference by adding reference in
        # 1st frame
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        ref = traj[5]
        act = adict['rmsd']
        mask = ":3-18@CA"
        dslist = act(mask, [[ref, ], traj], traj.top)
        saved_rmsd = traj.calc_rmsd(5, mask)
        # exclude 1st value for ref (=0.0)
        aa_eq(saved_rmsd, dslist[0].tolist()[1:])

if __name__ == "__main__":
    unittest.main()
