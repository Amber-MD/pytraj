from __future__ import absolute_import
from .Topology import Topology
from .DataSetList import DataSetList
from .DataFileList import DataFileList
from .ActionList import ActionList
from .AnalysisList import AnalysisList

class State:
    def __init__(self):
        self.actlist = ActionList()
        self.anllist = AnalysisList()
        self.dslist = DataSetList()
        self.dflist = DataFileList()
        self.trajlist = []
        self.toplist = []
        self.reflist = []

    def add_action(self, *args, **kwd):
        self.actlist.add_action(*args, **kwd)

    def add_reference(self, ref):
        self.reflist.append(ref)

    def add_top(self, top):
        self.toplist.append(top)

    def add_trajin(self, traj):
        if isinstance(traj, (list, tuple)):
            self.trajlist += traj
        else:
            self.trajlist.append(traj)

    def add_trajout(self, traj):
        pass

    def process(self):
        pass

    def run():
        pass

    def do_actions(self):
        pass

    def is_empty(self):
        # TODO : check
        condition = self.actlist.is_empty() and self.dslist.is_empty()
        condition += self.dflist.is_empty()

    def n_actions(self):
        return self.actlist.n_actions()

    def n_analyses(self):
        return self.anallist.n_analyses()

