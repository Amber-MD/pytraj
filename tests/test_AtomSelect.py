import unittest
import numpy as np
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.AtomSelect import AtomSelect

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        s = AtomSelect(traj=traj, top=traj.top, frameidx=5)

        # make sure to use `property` properly
        assert s.frameidx == 5
        s.frameidx = 1
        print(s.selected_frame[:1])
        assert s.frameidx == s._frameidx == 1
        s.frameidx = 1
        assert s.frameidx == s._frameidx == 1
        print(s.selected_frame.coords[:10])
        print(traj[1].coords[:10])
        assert_almost_equal(s.selected_frame.coords, traj[1].coords)
        s.frameidx = 9

        for i in range(9):
            s.frameidx = i
            assert_almost_equal(s.selected_frame.coords, traj[i].coords)

    def test_1(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        s = AtomSelect(traj=traj, top=traj.top)
        print(s.get_indices("@CA"))
        s.frameidx = 9
        print(s.get_indices("@H=").__len__())

    def test_2(self):
        print(dir(mdio))
        top, traj = mdio.loadpdb_rcsb("2kxc")
        farray = FrameArray()
        farray.top = top
        s = AtomSelect(traj=traj, top=top)

        #print s.get_indices("@CA")
        #mdio.writetraj(filename="./output/test_2kxc.pdb", top=top, traj=traj, indices=(0,), overwrite=True)
        rmsd, mat, v1, v2 = traj[0].rmsd(traj[1], get_mvv=True, top=traj.top, mask="@CA")
        print('rmsd = ', rmsd)
        rmsd, mat, v1, v2 = traj[0].rmsd(traj[1], get_mvv=True)
        print('rmsd = ', rmsd)
        print(mat, v1, v2)
        print(traj[0, 0])
        farray = traj[:]
        farray[0].rotate(mat)
        farray[0].translate(v1)
        print(farray[0, 0])

if __name__ == "__main__":
    unittest.main()
