import unittest
import os
import numpy as np
from pytraj.base import *
from pytraj.actions.Action_Molsurf import Action_Molsurf
from test_API.TestAPI import create_state, do_calculation

# not tested
#dir = os.environ["PYCPPTRAJ_HOME"] + "/tests/Cpptraj_test/"
#topdir = dir + "/tz2.parm7"
#crddir = dir + "/tz2.nc"
#state = create_state(top=topdir, trajin=crddir, ref=None)
#state.set_no_progress()
#surf = do_calculation(action=Action_Molsurf(), input="molsurf", state=state)
#
#print surf[:10]
#cppout = np.loadtxt(dir + "/Test_Molsurf/msurf.dat.save", skiprows=1).transpose()[1]
#print cppout[:10]
#np.testing.assert_almost_equal(surf[:10], cppout, decimal=3)
#print "Kool"
