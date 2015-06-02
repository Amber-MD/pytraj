from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having
import pytraj.common_actions as pyca

"""
try not to get segmentation fault error (due to whatever freaking reason)
"""
traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

class Test(unittest.TestCase):
    def test_0(self):
        print ("iter")
        from pytraj._shared_methods import _frame_iter_master
        it = _frame_iter_master(traj)

        for idx, frame in enumerate(it):
            pass
            #print ("segmentation faul if uncommenting #traj")
            #traj[idx]

        fa = traj[:]
        for idx, frame in enumerate(fa):
            fa[idx]

    def test_1(self):
        print ("calling search_hbonds several times")
        import pytraj.common_actions as pyca
        pyca.search_hbonds(traj)
        pyca.search_hbonds(traj, 'series')
        pyca.search_hbonds(traj, 'series, nointramol')

    def test_2(self):
        print ("DataSetList lifetime")
        #d = pyca.search_hbonds(traj)
        # FIXME: segmentation fault
        d = pyca.search_hbonds(traj).groupby("SER")
        d2 = pyca.search_hbonds(traj).groupby("SER").to_ndarray()
        print (d.size)
        print (d.keys())
        print (d)
        print (d2)

if __name__ == "__main__":
    unittest.main()
