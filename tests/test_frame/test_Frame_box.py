import unittest
import pytraj as pt
from utils import fn
import numpy as np


def test():
    traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
    frame0 = traj[0]
    assert frame0.has_box() == False
    frame0.box
    assert frame0.box.type == 'no_shape'
