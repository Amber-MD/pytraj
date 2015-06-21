from __future__ import print_function, absolute_import

from . DataSetList import DataSetList

class DSL(list):
    def __init__(self, ds=None):
        if ds:
            for d0 in ds:
                self.append(d0)
