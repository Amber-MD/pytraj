from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestClustering(unittest.TestCase):

    def test_kmeans(self):
        from pytraj import cluster

        kmeans = cluster.kmeans
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")

        cm = 'cluster kmeans clusters 10 @CA rms randompoint kseed 2 {0} noinfo'

        for sieve in [2, 3, 4]:
            sieve_str = "sieve {0} sieveseed 12".format(sieve)
            command = cm.format(sieve_str)
            state = pt.load_cpptraj_state(command, traj)
            state.run()

            data = kmeans(traj,
                          n_clusters=10,
                          kseed=2,
                          random_point=True,
                          metric='rms',
                          mask='@CA',
                          options=sieve_str)
            aa_eq(state.data[-1], data)


if __name__ == "__main__":
    unittest.main()
