#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.compat import PY3


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
        from pytraj.parallel.pjob import PJob
        traj = pt.load_sample_data('tz2')
        tasklist = []
        tasklist.append((pt.radgyr, traj))
        tasklist.append((pt.molsurf, traj, '@CA'))

        # perform each action on each CPUs (total 2 CPUs)
        pjob = PJob(tasklist)

        if PY3:
            data = pjob.compute()
            aa_eq(pt.radgyr(traj), data[0][1])
            aa_eq(pt.molsurf(traj, '@CA'), data[1][1])


if __name__ == "__main__":
    unittest.main()
    # nosetests --with-coverage --cover-package pytraj -vs .
    # nosetests -vs --processes 6 --process-timeout 200 .
    # nosetests -vs --processes 6 --process-timeout 200 --with-coverage --cover-package pytraj .
