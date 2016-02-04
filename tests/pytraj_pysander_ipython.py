from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
import numpy as np

try:
    import sander
    from parmed.amber.readparm import AmberParm
    has_sander_and_parmed = True
except:
    has_sander_and_parmed = False

if has_sander_and_parmed:
    traj_fn = "./data/Tc5b.x"
    top_fn = "./data/Tc5b.top"
    traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
    frame0 = traj[0]
    frame1 = traj[1]
    parm = AmberParm(top_fn)
    inp = sander.gas_input(8)
    parm.load_coordinates(frame0.coords)
else:
    print("require both sander and parmed. Skip test")
