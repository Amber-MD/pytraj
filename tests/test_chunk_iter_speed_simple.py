from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq, eq_coords
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        fa = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f_list = []

        for chunk in fa.chunk_iter(chunksize=2):
            f_list.append(chunk.average())

        print(f_list)
        f0 = fa[:2].average()
        aa_eq(f_list[0].xyz, f0.xyz)

        f0 = fa[2:4].average()
        aa_eq(f_list[1].xyz, f0.xyz)

if __name__ == "__main__":
    unittest.main()
