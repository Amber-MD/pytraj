from __future__ import print_function
from pytraj.actions import allactions

ADICT = {}

for key in allactions.__dict__.keys():
    if "Action_" in key:
        act = key.split('Action_')[1]
        #print (act)
        # create dict of action objects
        ADICT[act.lower()] = allactions.__dict__[key]

# make another dict to convert something like 'MolSurf' to 'molsurf'

class ActionDict:
    def __init__(self):
        self.adict = ADICT
        self.action_holder = None

    def __getitem__(self, key):
        # return Action object
        # why do we need action_holder?
        # should we use dict in command.cpp in cpptraj for mapping keyword
        # ('Action_DSSP' --> secstruct)
        self.action_holder = self.adict[key]()
        return self.action_holder

    def __del__(self):
        del self.action_holder

    def keys(self):
        return sorted(self.adict.keys())
