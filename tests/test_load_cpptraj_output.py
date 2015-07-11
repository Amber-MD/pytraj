from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


txt = """
parm ./data/Tc5b.top
trajin ./data/md1_prod.Tc5b.x
rms first"""


class Test(unittest.TestCase):

    def test_0(self):
        #pt.set_cpptraj_verbose()
        cout = pt.datafiles.load_cpptraj_output(txt)
        print(cout)

        cout = pt.datafiles.load_cpptraj_output(txt, with_state=True)
        print(cout)

if __name__ == "__main__":
    unittest.main()
