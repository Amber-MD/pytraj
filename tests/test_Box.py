import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.Box import Box
from array import array as pyarray

class TestBox(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        frame0.box_crd()
        print(frame0.get_box())
        frame0.boxview[:] = pyarray('d', [0.0, 1.0, 2.0, 3.0, 4.0, 6.])
        print(frame0.get_box())
        print(frame0.get_box().btype)
        frame0.set_nobox()
        print(frame0.get_box())

    def test_help(self):
        Box.help()

    def test_1(self):
        box = Box()
        box.set_trunc_oct()
        print(box)
        print(box.btype)
        box.set_nobox()
        print(box)
        print(box.btype)

if __name__ == "__main__":
    unittest.main()
