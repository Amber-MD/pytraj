from __future__ import print_function
import unittest

from pytraj import io as mdio
from pytraj.utils.check_and_assert import are_instance
from pytraj import *

from utils import fn


def test():
    traj = mdio.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
    assert are_instance([traj, traj], TrajectoryIterator) == True
    assert are_instance([traj, ""], TrajectoryIterator) == False
    assert are_instance(["my comment", ""], str) == True
