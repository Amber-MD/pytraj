import unittest
from pytraj import adict
import pytraj as pt
from utils import fn
from pytraj.datasets.c_datasetlist import DatasetList as CpptrajDatasetList


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))

        for i, frame in enumerate(traj.iterframe(1, 6, 2)):
            pass

        assert i == 2

        for i, frame in enumerate(traj.iterframe(1, 5, 1)):
            pass

        assert i == 3

        for i, frame in enumerate(traj.iterframe(stop=8)):
            pass

        assert i == 7

        for i, frame in enumerate(traj.iterframe(start=7, stop=8)):
            pass

        assert i == 0

    def test_1(self):
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        act = adict['distance']
        dslist = CpptrajDatasetList()
        act.read_input(":2@CA :10@CA", traj.top, dslist=dslist)
        act.setup(traj.top)

        for frame in traj.iterframe(stop=5):
            act.compute(frame)

        dslist = CpptrajDatasetList()
        act2 = adict['distance']
        act2.read_input(":2@CA :10@CA", traj.top, dslist=dslist)
        act2.setup(traj.top)
        act2.compute(traj.iterframe(stop=5))
        assert act2.n_frames == 5


if __name__ == "__main__":
    unittest.main()
