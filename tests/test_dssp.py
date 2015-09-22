import unittest
import numpy as np
import pytraj as pt
from pytraj.testing import aa_eq



class TestDSSP(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        data = pt.dssp(traj, "*", dtype='cpptraj_dataset')
        data_int = np.array([d0.values for d0 in data if d0.dtype == 'integer'], dtype='i4')
        print(data_int)

        # load cpptraj output
        cpp_data = np.loadtxt("./data/dssp.Tc5b.dat", skiprows=1).transpose()[1:]
        print(cpp_data)
        aa_eq(data_int.flatten(), cpp_data.flatten())


if __name__ == "__main__":
    unittest.main()
