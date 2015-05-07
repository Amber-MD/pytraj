from __future__ import absolute_import
import os
from .. TrajectoryIterator import TrajectoryIterator

def load_sample_data(data_name=None):
    """
    Return TrajectoryIterator instance for Ala3 or tz2 data

    Paramters
    ---------
    data_name : str, 'ala3' (defaul) | 'tz2'

    Notes
    -----
    tz2 dataset : $AMBERHOME/AmberTools/test/cpptraj/
        explicit water, ortho box
    """
    mydir = os.path.dirname(os.path.abspath(__file__))
    if data_name is None or data_name.lower() == 'ala3':
        crd = os.path.join(mydir, "Ala3", "Ala3.crd")
        top = os.path.join(mydir, "Ala3", "Ala3.top")
    elif data_name.lower() == 'tz2':
        crd = os.path.join(mydir, "tz2", "tz2.ortho.nc")
        top = os.path.join(mydir, "tz2", "tz2.ortho.parm7")
    return TrajectoryIterator(crd, top)
