import pytraj as pt
from pytraj.testing import aa_eq
import numpy as np


def test_hausdorff(tmpdir):
    # with tmpdir.as_cwd():
    matrix = np.array(
        [[1.0000, 2.2361, 3.0000, 4.1231], [2.2361, 1.0000, 2.2361, 5.0000],
         [3.0000, 2.2361, 1.0000, 4.1231], [2.2361, 3.0000, 2.2361, 3.0000]])
    aa_eq(pt.hausdorff(matrix), [[3.][3.][2.23609996]])
