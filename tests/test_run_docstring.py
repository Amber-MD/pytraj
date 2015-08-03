from __future__ import print_function
import unittest
from pytraj.testing import run_docstring
import pytraj.common_actions as pyca
from pytraj.base import *
import pytraj as pt


def silly_doc_func():
    """
    >>> print (traj)
    >>> print ("work nicely")
    """
    pass


class Test(unittest.TestCase):

    def test_0(self):
        print("silly_doc_func")
        run_docstring(silly_doc_func)

        print("_frame_iter_master")
        from pytraj._shared_methods import _frame_iter_master as fi
        run_docstring(fi)

        print("matrix_analysis")
        from pytraj import matrix_analysis as ma
        func_names = ma.default_key_dict.keys()
        for key in func_names:
            run_docstring(ma.__dict__[key])

        print("dihedral_analysis")
        from pytraj import dihedral_analysis as da
        func_names = da.supported_dihedral_types
        for key in func_names:
            run_docstring(da.__dict__['calc_' + key])

        run_docstring(pt.multidihedral)

        print("pyca.calc_rmsd")
        run_docstring(pyca.calc_rmsd)

        print("Trajectory.calc_rmsd")
        from pytraj import Trajectory
        run_docstring(Trajectory.calc_rmsd)

        print("Topology.select")
        run_docstring(Topology.select)

        print("frame __getitem__")
        from pytraj import Frame
        run_docstring(Frame.__getitem__)

        print("split_range")
        from pytraj.misc import split_range
        run_docstring(split_range)

        print("pytraj.tools.grep")
        from pytraj.tools import grep
        run_docstring(grep)

        print ("load_ParmEd")
        run_docstring(pt.load_ParmEd)

        print("load Topology")
        run_docstring(pt.load_topology)

        print("clustering_dataset")
        from pytraj import clustering_dataset
        pt.set_cpptraj_verbose()
        run_docstring(clustering_dataset)
        pt.set_cpptraj_verbose(False)

if __name__ == "__main__":
    unittest.main()
