from __future__ import absolute_import
from pytraj.datasets.DataSetList import DataSetList as DSL
from pytraj.externals._json import to_json, read_json
from pytraj.externals._pickle import to_pickle, read_pickle


class DataSetList(DSL):
    def from_pickle(self, filename):
        from pytraj.datasets import cast_dataset
        ddict = read_pickle(filename)
        ordered_keys = ddict['ordered_keys']

        for legend in ordered_keys:
            d = ddict[legend]
            values = d['values']
            self.add_set(d['dtype'], d['name'])
            last = self[-1]
            last.set_name_aspect_index_ensemble_num(d['aspect'], d['name'], d['idx'], 0)
            last.set_legend(legend)
            last.resize(len(values))
            last.values[:] = values

    def to_pickle(self, filename, use_numpy=True):
        to_pickle(self._to_full_dict(use_numpy), filename)

    def _to_full_dict(self, use_numpy=True):
        """
        """
        ddict  = {}
        ddict['ordered_keys'] = []
        for d in self:
            ddict['ordered_keys'].append(d.legend)
            ddict[d.legend] = {}
            _d = ddict[d.legend]
            if use_numpy:
                _d['values'] = d.values
            else:
                _d['values'] = list(d.values)
            _d['name'] = d.name
            _d['dtype'] = d.dtype
            _d['aspect'] = d.aspect
            _d['idx'] = d.idx
        return ddict

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
