from __future__ import print_function
import unittest
from pytraj.testing import run_docstring
import pytraj.common_actions as pyca

def silly_doc_func():
    """
    >>> print (traj)
    >>> print ("work nicely")
    """
    pass

class Test(unittest.TestCase):
    def test_0(self):
        print ("silly_doc_func")
        run_docstring(silly_doc_func)

        print ("_frame_iter_master")
        from pytraj._shared_methods import _frame_iter_master as fi
        run_docstring(fi)

        print ("matrix_analysis")
        from pytraj import matrix_analysis as ma
        func_names = ma.default_key_dict.keys()
        for key in func_names:
            run_docstring(ma.__dict__[key])

        print ("pyca.calc_rmsd")
        run_docstring(pyca.calc_rmsd)

        print ("Trajectory.calc_rmsd")
        from pytraj import Trajectory
        run_docstring(Trajectory.calc_rmsd)

if __name__ == "__main__":
    unittest.main()
