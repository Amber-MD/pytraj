from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq, a_isinstance
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    @test_if_having("numpy")
    def test_0(self):
        import numpy as np
        from array import array as pyarray
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        alist = traj.calc_molsurf().tolist()
        anp = traj.calc_molsurf().to_ndarray()
        a_pyarray = traj.calc_molsurf().to_pyarray()
        aa_eq(alist, anp)
        aa_eq(alist, a_pyarray)

        alist2 = pyca.calc_distance(traj, "@2 @4", dtype='list')
        anp2 = pyca.calc_distance(traj, "@2 @4", dtype='ndarray')
        a_pyarray2 = pyca.calc_distance(traj, "@2 @4", dtype='pyarray')
        a_isinstance(alist2, list)
        a_isinstance(anp2, np.ndarray)
        a_isinstance(a_pyarray2, pyarray)

        aa_eq(alist2, anp2)
        aa_eq(alist2, a_pyarray2)

        # test hist
        d0 = traj.calc_molsurf()
        print (d0.hist(bins=3, range=[d0.min(), d0.max()]))
        print (d0.to_ndarray())

if __name__ == "__main__":
    unittest.main()
