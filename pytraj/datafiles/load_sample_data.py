from __future__ import absolute_import
import os
from ..TrajectoryIterator import TrajectoryIterator


def load_sample_data(data_name=None):
    """
    Return TrajectoryIterator instance for Ala3 or tz2 data

    Paramters
    ---------
    data_name : str, {'ala3', 'tz2', 'rna'}, default 'ala3'

    Notes
    -----
    tz2 dataset : $AMBERHOME/AmberTools/test/cpptraj/
        explicit water, ortho box
    """
    data_dict = {'ala3': ["Ala3/Ala3.crd", "Ala3/Ala3.top"],
                 'tz2' : ["tz2/tz2.ortho.nc", "tz2/tz2.ortho.parm7"],
                 'rna' : ["rna.pdb", "rna.pdb"]
                 }

    mydir = os.path.dirname(os.path.abspath(__file__))
    if data_name is None:
        data_name = 'ala3'
    crd = os.path.join(mydir, data_dict[data_name][0])
    top = os.path.join(mydir, data_dict[data_name][1])
    return TrajectoryIterator(crd, top)

def load_rna():
    '''return pytraj.TrajectoryIterator for an RNA trajectory with 3 frames
    '''
    return load_sample_data('rna')
