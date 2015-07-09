
import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):

    def test_0(self):
        fname = './data/Tc5b.top'
        top = Topology(fname)
        top.summary()

    def test_1(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        print(traj)
        print(traj.top)
        print(traj.n_frames)
        print(traj[0, 0, :])
        traj[:4].save('./output/test_python3.x', overwrite=True)

    def test_2(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        top = traj.top
        print(top("@CA"))
        print(top["@CA"][:2])
        atom = top["@CA"][0]
        print(atom.name)
        print(type(atom.name))
        print(atom.name == 'CA  ')
        top.atom_info('@CA')

if __name__ == "__main__":
    unittest.main()
