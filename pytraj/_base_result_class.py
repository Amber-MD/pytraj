
class BaseAnalysisResult(object):
    def __init__(self, dslist):
        self.dslist = dslist

    def to_ndarray(self):
        return self.dslist.to_ndarray()
