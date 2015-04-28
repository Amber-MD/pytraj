from __future__ import print_function
import unittest
from pytraj import Topology

class Test(unittest.TestCase):
    def test_0(self):
        top = Topology("./data/tz2.ortho.parm7")
        box = top.box
        print (box)
        from pytraj.misc import get_atts

        for att in get_atts(box):
            print (getattr(box, att))

if __name__ == "__main__":
    unittest.main()
