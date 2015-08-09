from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        from pytraj.testing import load_sample_data
        from pytraj.utils import _import_pandas
        traj = load_sample_data()

        _, pd = _import_pandas()
        if not pd:

            def test():
                import pandas

            self.assertRaises(ImportError, lambda: test())
        else:
            df = traj.top.to_dataframe()
            print(df.__str__())


if __name__ == "__main__":
    unittest.main()
