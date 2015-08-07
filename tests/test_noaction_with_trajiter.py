from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        self.assertRaises(ValueError, lambda: pyca.translate(traj, '@1 1.0'))
        self.assertRaises(ValueError, lambda: pyca.rotate(traj, '@1 1.0'))
        self.assertRaises(ValueError, lambda: pyca.autoimage(traj, '@1 1.0'))
        self.assertRaises(ValueError, lambda: pyca.scale(traj, '@1 1.0'))

        traj2 = traj.to_mutable_trajectory()
        pyca.translate(traj2, '@1 1.0')
        pyca.rotate(traj2, 'x 1.0')
        pyca.autoimage(traj2)
        pyca.scale(traj2, 'x 0.5')


if __name__ == "__main__":
    unittest.main()
