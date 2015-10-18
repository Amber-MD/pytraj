from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
import pytraj.common_actions as pyca

try:
    import sander
    has_sander = True
except ImportError:
    has_sander = False

@unittest.skipIf(not has_sander, 'skip if not having sander')
class TestSanderPmap(unittest.TestCase):
    def test_sander_pmap(self):
        traj = pt.iterload('./data/md1_prod.Tc5b.x', './data/Tc5b.top')
        fname = traj.top.filename
        print(pt.energy_decomposition(traj, parm=fname)['dihedral'])
        print([x[1]['dihedral'] for x in pt.pmap(n_cores=4,
                                                 func=pt.energy_decomposition,
                                                 traj=traj,
                                                 parm=fname)])


if __name__ == "__main__":
    unittest.main()
