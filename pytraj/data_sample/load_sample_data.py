from __future__ import absolute_import
import os
from .. import io as mdio

def load_sample_data():
    """
    Return FrameArray instance for Ala3 data
    """
    mydir = os.path.dirname(os.path.abspath(__file__))
    # FIXME: update this if want to support Window
    return mdio.load(mydir+"/Ala3/Ala3.crd", mydir+"/Ala3/Ala3.top")
