import unittest
from array import array
from pytraj.base import *
import numpy as np

class TestCHARMM(unittest.TestCase):
    def test_0(self):
        top = Topology("./data/ala3.psf")
        print(top)
        print(top.n_atoms)
        print(dir(top))
        top.summary()
        reslit = top.residuelist
        print(reslit[0])
        atm = AtomMask("@CA")
        print(atm.n_selected)
        print(atm)
        print(dir(atm))
        top.set_integer_mask(atm)
        print(atm.n_selected)
        print(atm.atoms_in_char_mask())
        atm.invert_mask()
        atm.mask_info()
        atm.invert_mask()
        atm.mask_info()

        atm.invert_mask()
        frame = Frame(atm.n_selected)
        frame[:10] = np.asarray(array('d', list(range(30)))).reshape(10, 3)
        print(frame[0])
        print(dir(frame))
        frame.zero_coords()
        print(frame[0])
        print(atm.is_char_mask())

if __name__ == "__main__":
    unittest.main()
