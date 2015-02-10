from pytraj.Topology import Topology
from pytraj.DataSetList import DataSetList
from pytraj.DataFileList import DataFileList
from pytraj.ActionList import ActionList

class State:
    def __init__(self):
        self._actlist = ActionList()
        self._dslist = DataSetList()
        self._dflist = DataFileList()
        self.top = Topology()

    def add_action(self, *args, **kwd):
        self._actlist.add_action(*args, **kwd)

    def process(self):
        pass

    def do_actions(self):
        pass

    def is_empty(self):
        # TODO : check
        condition = self._actlist.is_empty() and self._dslist.is_empty()
        condition += self._dflist.is_empty()

    def n_actions(self):
        return self._actlist.n_actions()
