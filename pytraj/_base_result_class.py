from pytraj.datasetlist import _groupby
from pytraj.datasetlist import DatasetList


class BaseAnalysisResult(object):
    def __init__(self, dslist=None):
        if dslist is not None:
            self.dslist = dslist
        else:
            self.dslist = DatasetList()

    def to_ndarray(self):
        return self.dslist.to_ndarray()

    def to_dict(self):
        return self.dslist.to_dict()

    def grep(self, key):
        return self.dslist.grep(key)

    def __getitem__(self, idx):
        return self.__class__(self.dslist[idx])

    def __iter__(self):
        return self.dslist.__iter__()

    def groupby(self, key):
        return _groupby(self, key)

    def append(self, value):
        self.dslist.append(value)

    @property
    def values(self):
        return self.dslist.values
