from __future__ import absolute_import
from pkg_resources import resource_filename
from .. import io as mdio

def load_sample_data():
    """
    Return FrameArray instance for Ala3 data
    """
    mydir = resource_filename('pytraj', 'data_sample/')
    return mdio.load(mydir+"Ala3/Ala3.crd", mydir+"./Ala3/Ala3.top")
