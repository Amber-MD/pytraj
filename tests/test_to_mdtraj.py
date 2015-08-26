from __future__ import print_function
import unittest; import pytraj as pt
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    @test_if_having("mdtraj")
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        traj2 = pt.load_mdtraj(pt.to_mdtraj(traj), autoconvert=False)
        print(traj[0, 0], traj2[0, 0])
        aa_eq(traj.xyz, traj2.xyz)


if __name__ == "__main__":
    unittest.main()
