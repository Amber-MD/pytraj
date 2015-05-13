from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from pytraj.utils import Timer
from pytraj.externals.six.moves import range

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        f0 = traj[0]
        f0.strip_atoms(mask='!@CA', top=traj.top)
        atm0 = traj.top("!@CA")
        atm0.invert_mask()
        NATOM = atm0.n_atoms
        assert f0.n_atoms == NATOM
        print (f0)

        fa0 = traj[:]
        fa0.strip_atoms('!@CA')
        assert fa0[0].n_atoms == NATOM
        fa1 = traj[:]
        fa1._fast_strip_atoms('!@CA')
        assert fa1[0].n_atoms == NATOM

        aa_eq(fa0.xyz, fa1.xyz)

        #fa3 = traj[:]
        fa3 = traj._fast_slice(slice(2, 100, 10))
        print (traj)
        #print (fa3)
        #for i in range(20):
        #    print (fa3)
        #    fa3 += fa3.copy()

        #fa4 = fa3.copy()
        #print (fa3, fa4)

        #@Timer()
        #def normal_strip():
        #    fa3.strip_atoms("!@H,C")

        #@Timer()
        #def openmp_strip():
        #    fa4._fast_strip_atoms("!@H,C")

        #normal_strip()
        #openmp_strip()
        #aa_eq(fa3.xyz, fa4.xyz)
        #atm = traj.top("!@H,C")
        #atm.invert_mask()
        #assert fa3.n_atoms == atm.n_atoms == fa4.n_atoms

if __name__ == "__main__":
    unittest.main()
