
class BaseAnalysisResult(object):

    def __init__(self, dslist):
        self.dslist = dslist

    def to_ndarray(self):
        return self.dslist.to_ndarray()

    def to_dict(self):
        return self.dslist.to_dict()

    def grep(self, key):
        return self.__class__(self.dslist.grep(key))

    def __getitem__(self, idx):
        return self.__class__(self.dslist[idx])

    def __iter__(self):
        return self.dslist.__iter__()
