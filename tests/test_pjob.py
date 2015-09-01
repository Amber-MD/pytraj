from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.tools import PY2


class Test(unittest.TestCase):
    @unittest.skipIf(PY2, 'only work python3')
    def test_0(self):
        from pytraj.parallel import PJob
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        job = PJob([(pt.radgyr, traj), (pt.molsurf, traj)])
        results = job.compute()
        #print(results)

        aa_eq(pt.radgyr(traj), results[0][1])
        aa_eq(pt.molsurf(traj), results[1][1])


if __name__ == "__main__":
    unittest.main()
