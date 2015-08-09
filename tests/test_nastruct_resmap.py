from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        f0 = "./data/nastruct/rna_a_form_dup.pdb"
        f1 = "./data/nastruct/rna_a_form_dup_F2.pdb"
        traj0 = pt.iterload(f0)
        traj1 = pt.iterload(f1)
        print(traj0)

        a0 = pt.nastruct(traj0)
        a1 = pt.nastruct(traj1, resmap='AF2:A')
        aa_eq(a0.values, a1.values)
        print(len(a0.values))

        cout = pt.datafiles.load_cpptraj_output("""
        parm ./data/nastruct/rna_a_form_dup.pdb
        trajin ./data/nastruct/rna_a_form_dup.pdb
        nastruct naout na.out
        """,
                                                dtype='cpptraj_dataset',
                                                with_traj=True)
        print(cout)
        aa_eq(traj0.xyz, cout[0].xyz)

        expected_result = cout[1]
        print(expected_result.to_ndarray())


if __name__ == "__main__":
    unittest.main()
