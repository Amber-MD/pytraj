from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        from pytraj import cluster
        pt.set_cpptraj_verbose()

        kmeans = cluster.kmeans
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        dset = kmeans(traj,
                      n_clusters=10,
                      kseed=2,
                      random_point=True,
                      distance_metric='rms',
                      mask='@CA')


if __name__ == "__main__":
    unittest.main()
