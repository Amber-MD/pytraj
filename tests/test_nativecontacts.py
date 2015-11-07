#!/usr/bin/env python
from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestNativeContacts(unittest.TestCase):
    def test_nativecontacts(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        dslist = pt.native_contacts(traj, top=traj.top)
        cpp = np.loadtxt('data/tc5b.native_contacts.dat',
                         skiprows=1,
                         usecols=(1, 2)).T
        aa_eq(dslist.values, cpp)


if __name__ == "__main__":
    unittest.main()
