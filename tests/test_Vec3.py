from pytraj import Vec3
#from pytraj import DistRoutines as dist
import numpy as np
import unittest

class Test(unittest.TestCase):
    def test_1(self):
        v1 = Vec3(0., 0., 0.)
        print(v1[:])
        v1[0] = 200.
        v1[:] = [100, 300, 400]
        assert v1.tolist() == [100, 300, 400]
        assert isinstance(v1.tolist(), list) == True
        print(np.asarray(v1[:]))
        print(v1)

        from pytraj.math import Matrix_3x3
        mat = Matrix_3x3(list(range(9)))
        print(mat*v1)

if __name__ == "__main__":
    unittest.main()
