from __future__ import absolute_import, print_function
from pytraj.base import *
from pytraj.core import *
from pytraj.math import *
from pytraj.externals import *
from pytraj.c_traj import *
from pytraj.hbond_analysis import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from .datafiles.load_samples import load_sample_data
from pytraj.misc import info
from pytraj.c_options import set_world_silent
from pytraj import trajectory

from pytraj import *
from pytraj.datasets import *
from pytraj.all_actions import *
from pytraj.c_action import c_action
from pytraj.c_analysis import c_analysis


def run_tests():
    print("try to load sample data")
    traj = load_sample_data()
    traj = load_sample_data('tz2')

    print("try to make all action objects")
    failed_list = ['createreservoir', ]
    DatasetList()
    print("try to make all analysis objects")
    from pytraj import analdict
    failed_list = []

    for key in analdict.keys():
        if key not in failed_list:
            analdict[key]

    print("try to make all dataset stuff")
    DatasetDouble()
    DatasetFloat()
    DatasetInteger()
    DatasetString()
    DatasetMatrixDouble()
    DatasetGridFloat()
    DatasetMatrixFloat()
    DatasetVector()
    DatasetMatrix3x3()
    DatasetCoords()
    DatasetCoordsRef()
    DatasetCoordsCRD()

    print("try to make structure-related objects")
    Topology()
    Molecule()
    Residue()
    Atom()
    Box()
    Frame()

    print("try to create Trajectory-like objects")
    TrajectoryIterator()
    trajectory.Trajectory()

    print("other stuff. throw all tests don't belong anywhere else here")
    from pytraj import c_dict
    from pytraj.misc import get_atts
    keys = get_atts(c_dict)
    cdict = c_dict.__dict__

    for key in keys:
        if isinstance(cdict[key], dict):
            assert cdict[key].keys() is not None

    # other objects
    CpptrajState()

    print("OK")


if __name__ == '__main__':
    run_tests()
