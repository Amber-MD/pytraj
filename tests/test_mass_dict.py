from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):

    def test_0(self):
        from pytraj.core import mass_atomic_number_dict, mass_element_dict
        top = mdio.load_topology("./data/tz2.parm7")
        mass_list = []

        for atom in top:
            mass_list.append(mass_atomic_number_dict[atom.atomic_number])
        aa_eq(mass_list, top.mass, 2)


if __name__ == "__main__":
    unittest.main()
