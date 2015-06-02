from __future__ import print_function
import unittest
import sys
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having
import pytraj.common_actions as pyca


class TestReferentCounting(unittest.TestCase):
    def test_pyca(self):
        # pyca
        traj = mdio.iterload('./data/md1_prod.Tc5b.x', './data/Tc5b.top')
        
        # pyca
        d = pyca.search_hbonds(traj)
        
        assert sys.getrefcount(d) == 2
        sys.getrefcount(d) == 2
        junk = d.groupby('SER'); del junk
        
        assert sys.getrefcount(d) == 2
        sys.getrefcount(d)
        
        pyca.search_hbonds(traj).groupby('SER')
        pyca.search_hbonds(traj).groupby('SER').groupby("SER")
        e = pyca.search_hbonds(traj).groupby('SER').groupby("SER").groupby("").groupby("")
        # make sure getting not segmentation fault
        print (e.size)
        print (e[0])

    def test_traj_search_hbonds(self):
        traj = mdio.iterload('./data/md1_prod.Tc5b.x', './data/Tc5b.top')
        # traj.search_hbonds
        d = traj.search_hbonds()
        
        assert sys.getrefcount(d) == 2
        sys.getrefcount(d) == 2
        junk = d.groupby('SER'); del junk
        
        assert sys.getrefcount(d) == 2
        sys.getrefcount(d)
        
        traj.search_hbonds().groupby('SER')
        traj.search_hbonds().groupby('SER').groupby("SER")
        e = traj.search_hbonds().groupby('SER').groupby("SER").groupby("").groupby("")
        # make sure getting not segmentation fault
        print (e.size)
        print (e[0])

if __name__ == '__main__':
    unittest.main()
