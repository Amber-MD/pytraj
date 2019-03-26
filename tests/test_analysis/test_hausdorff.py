import pytraj as pt
import numpy as np
from numpy.testing import assert_almost_equal as aa_eq


def test_hausdorff(tmpdir):
    # with tmpdir.as_cwd():
    matrix = np.array(
        [[1.0000, 2.2361, 3.0000, 4.1231], [2.2361, 1.0000, 2.2361, 5.0000],
         [3.0000, 2.2361, 1.0000, 4.1231], [2.2361, 3.0000, 2.2361, 3.0000]])
    aa_eq(pt.hausdorff(matrix), [[3.], [3.], [2.23609996]])
