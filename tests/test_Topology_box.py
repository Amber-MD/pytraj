from __future__ import print_function
import unittest
from pytraj.testing import aa_eq
from pytraj import Topology


class Test(unittest.TestCase):
    def test_0(self):
        top = Topology("./data/tz2.ortho.parm7")
        box = top.box
        #print(box)
        from pytraj.misc import get_atts

        for att in get_atts(box):
            #print(getattr(box, att))

        top2 = Topology()
        box = top.box
        # list
        top2.box = box.tolist()
        aa_eq(top2.box.tolist(), top.box.tolist())
        assert top2.box.type == top.box.type

        # tuple
        top2.box = tuple(box.tolist())
        aa_eq(top2.box.tolist(), top.box.tolist())
        assert top2.box.type == top.box.type

        # memview
        top2.box = box[:]
        aa_eq(top2.box.tolist(), top.box.tolist())
        assert top2.box.type == top.box.type


if __name__ == "__main__":
    unittest.main()
