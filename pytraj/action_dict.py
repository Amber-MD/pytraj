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

    def __getitem__(self, key):
        # return Action object
        return self.adict[key]()

    def keys(self):
        return self.adict.keys()
