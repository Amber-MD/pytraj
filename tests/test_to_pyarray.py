from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists

class Test(unittest.TestCase):
    @test_if_having("numpy")
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        alist = traj.calc_molsurf().tolist()
        anp = traj.calc_molsurf().to_ndarray()
        a_pyarray = traj.calc_molsurf().to_pyarray()
        aa_eq(alist, anp)
        aa_eq(alist, a_pyarray)

if __name__ == "__main__":
    unittest.main()
