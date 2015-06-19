from __future__ import absolute_import, print_function
from pytraj.base import *
from pytraj.core import *
from pytraj.math import *
from pytraj.externals import *
from pytraj.plotting import *
from pytraj.trajs import *
from pytraj.hbonds import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from .data_sample.load_sample_data import load_sample_data
from pytraj import adict
from pytraj.misc import info
from pytraj._set_silent import set_world_silent
from pytraj import api

from pytraj import *
from pytraj.parms import *
from pytraj.datasets import *
from pytraj.common_actions import *
from pytraj.actions import CpptrajActions
from pytraj.analyses import CpptrajAnalyses

def run_tests():
    print ("try to load sample data")
    traj = load_sample_data()
    traj = load_sample_data('tz2')

    print ("try to make all action objects")
    from pytraj import adict
    # FIXME, TODO : make failed_list empty
    failed_list = ['createreservoir',]

    for key in adict.keys():
        if key not in failed_list:
            adict[key]

    DataSetList()
    print ("try to make all analysis objects")
    from pytraj import analdict
    failed_list = []

    for key in analdict.keys():
        if key not in failed_list:
            analdict[key]

    print ("try to make all dataset stuff")
    DatasetDouble()
    DatasetFloat()
    DatasetInteger()
    DatasetString()
    DatasetMatrixDouble()
    DatasetGridFloat()
    DataSet_MatrixFlt()
    DatasetVector()
    DataSet_Coords()
    DataSet_Coords_REF()
    DataSet_Coords_CRD()
    DataSet_Coords_TRJ()

    print ("try to make structure-related objects")
    Topology()
    Molecule()
    Residue()
    Atom()
    Box()
    Frame()

    print ("try to create Trajectory-like objects")
    Trajectory()
    TrajectoryIterator()
    TrajectoryREMDIterator.TrajectoryREMDIterator()
    TrajinList.TrajinList()
    api.Trajectory()

    print ("other stuff. throw all tests don't belong anywhere else here")
    from pytraj import cpptraj_dict
    from pytraj.misc import get_atts
    keys = get_atts(cpptraj_dict)
    cdict = cpptraj_dict.__dict__

    for key in keys:
        if isinstance(cdict[key], dict):
            assert cdict[key].keys() is not None

    # other objects
    CpptrajState()

    print ("OK")

if __name__ == '__main__':
    run_tests()
