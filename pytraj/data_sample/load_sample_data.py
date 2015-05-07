from __future__ import absolute_import
import os
from .. TrajectoryIterator import TrajectoryIterator

def load_sample_data():
    """
    Return TrajectoryIterator instance for Ala3 data
    """
    mydir = os.path.dirname(os.path.abspath(__file__))
    crd = os.path.join(mydir, "Ala3", "Ala3.crd")
    top = os.path.join(mydir, "Ala3", "Ala3.top")
    return TrajectoryIterator(crd, top)
