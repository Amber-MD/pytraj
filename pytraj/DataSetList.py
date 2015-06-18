from __future__ import absolute_import
from pytraj.datasets.DataSetList import DataSetList as DSL
from pytraj.externals._json import to_json


class DataSetList(DSL):
    def to_json(self, filename):
        to_json(self.to_dict(), filename)

    def to_dataframe(self, engine='pandas'):
        if engine == 'pandas':
            try:
                import pandas as pd
                return pd.DataFrame(self.to_dict(use_numpy=True))
            except ImportError:
                raise ImportError("must have pandas")
        else:
            raise NotImplementedError("currently support only pandas' DataFrame")

    def hist(self, plot=False):
        """
        Paramters
        ---------
        plot : bool, default False
            if False, return a dictionary of 2D numpy array
            if True, return a dictionary of matplotlib object
        """
        return dict(map(lambda x : (x.legend,  x.hist(plot=plot)), self))

    def count(self):
        from collections import Counter
        return dict((d0.legend, Counter(d0.values)) for d0 in self)

    def chunk_average(self, n_chunks):
        return dict((d0.legend, d0.chunk_average(n_chunks)) for d0 in self)

    def dtypes(self):
        return self.get_dtypes()

    def aspects(self):
        return self.get_aspects()

    def pipe(self, *funcs):
        """apply a series of functions to self's data
        """
        values = self.values
        for func in funcs:
            values = func(values)
        return values
