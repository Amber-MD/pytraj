from __future__ import print_function
import unittest
import pytraj as pt
import pytraj as pt
from pytraj.testing import run_docstring
import pytraj.common_actions as pyca
from pytraj.base import *
import pytraj as pt


def silly_doc_func():
    """
    """
    pass


class Test(unittest.TestCase):
    def test_0(self):
        for func in [pt.iterframe,
                     pt.iterchunk,
                     pt.watershell,
                     pt.distance,
                     pt.angle,
                     pt.dihedral,
                     pt.dssp,
                     pt.radgyr,
                     pt.molsurf,
                     pt.search_hbonds,
                     pt.closest,
                     pt.search_neighbors,
                     pt.pairwise_rmsd, ]:
            run_docstring(func)

        run_docstring(silly_doc_func)

        from pytraj._shared_methods import _frame_iter_master as fi
        run_docstring(fi)

        from pytraj import matrix_analysis as ma
        func_names = ma.default_key_dict.keys()
        for key in func_names:
            run_docstring(ma.__dict__[key])

        from pytraj import dihedral_analysis as da
        func_names = da.supported_dihedral_types
        for key in func_names:
            run_docstring(da.__dict__['calc_' + key])

        run_docstring(pt.multidihedral)

        run_docstring(pyca.calc_rmsd)

        from pytraj import Trajectory
        run_docstring(Topology.select)

        from pytraj import Frame
        run_docstring(Frame.__getitem__)

        from pytraj.misc import split_range
        run_docstring(split_range)

        from pytraj.tools import grep
        run_docstring(grep)

        run_docstring(pt.load_ParmEd)

        run_docstring(pt.load_topology)

        from pytraj import clustering_dataset
        run_docstring(clustering_dataset)
        run_docstring(pt.mindist)

        from pytraj.cluster import kmeans
        run_docstring(kmeans)


if __name__ == "__main__":
    unittest.main()
