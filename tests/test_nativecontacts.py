#!/usr/bin/env python
from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestNativeContacts(unittest.TestCase):

    def test_nativecontacts(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        dslist = pt.native_contacts(traj, top=traj.top)
        cpp = np.loadtxt('data/tc5b.native_contacts.dat',
                         skiprows=1,
                         usecols=(1, 2)).T
        aa_eq(dslist.values, cpp)

        # mask2
        cb_indices = pt.select('@CB', traj.top)
        dslist2 = pt.native_contacts(traj, mask='@CA', mask2='@CB')
        dslist3 = pt.native_contacts(traj, mask='@CA', mask2=cb_indices)
        aa_eq(dslist2.values, dslist3.values)


if __name__ == "__main__":
    unittest.main()
