from __future__ import absolute_import
from pytraj.datasets.DataSetList import DataSetList as DSL
from pytraj.externals._json import to_json


class DataSetList(DSL):
    def to_json(self, filename):
        to_json(self.to_dict(), filename)

    def to_dataframe(self):
        try:
            import pandas as pd
            return pd.DataFrame(self.to_dict(use_numpy=True))
        except ImportError:
            raise ImportError("must have pandas")

    def hist(self, plot=False):
        """
        Paramters
        ---------
        plot : bool, default False
            if False, return a dictionary of 2D numpy array
            if True, return a dictionary of matplotlib object
        """
        return dict(map(lambda x : (x.legend,  x.hist(plot=plot)), self))
