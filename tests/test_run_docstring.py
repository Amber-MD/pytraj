from __future__ import print_function
import unittest
from pytraj.testing import run_docstring

class Test(unittest.TestCase):
    def test_0(self):
        # _frame_iter_master
        from pytraj._shared_methods import _frame_iter_master as fi
        run_docstring(fi)

        # matrix_analysis
        from pytraj import matrix_analysis as ma
        func_names = ma.default_key_dict.keys()
        for key in func_names:
            run_docstring(ma.__dict__[key])

if __name__ == "__main__":
    unittest.main()
