from __future__ import absolute_import
from pytraj.datasets.DataSetList import DataSetList as DSL
from pytraj.externals._json import to_json

# subclass DataSetList so we can use in python level

#class DataSetList(object):
class DataSetList(DSL):
    def to_json(self, filename):
        to_json(self.to_dict(), filename)

    def to_dataframe(self):
        try:
            import pandas as pd
            return pd.DataFrame(self.to_dict(use_numpy=True))
        except ImportError:
            raise ImportError("must have pandas")
