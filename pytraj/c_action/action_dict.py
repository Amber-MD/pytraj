from __future__ import print_function
from . import c_action as allactions

ADICT = {}

for key in allactions.__dict__.keys():
    if "Action_" in key:
        act = key.split('Action_')[1]
        # create dict of action objects
        ADICT[act.lower()] = allactions.__dict__[key]

# add some commond words to ADICT
ADICT['surf_LCPO'] = allactions.__dict__['Action_Surf']
ADICT['surf_lcpo'] = allactions.__dict__['Action_Surf']
ADICT['secstruct'] = allactions.__dict__['Action_DSSP']
ADICT['rms'] = allactions.__dict__['Action_Rmsd']
ADICT['superpose'] = allactions.__dict__['Action_Rmsd']
ADICT['drmsd'] = allactions.__dict__['Action_DistRmsd']
ADICT["lipidorder"] = allactions.__dict__['Action_OrderParameter']
ADICT["rog"] = allactions.__dict__['Action_Radgyr']
ADICT["stfcdiffusion"] = allactions.__dict__['Action_STFC_Diffusion']
ADICT["symmrmsd"] = allactions.__dict__['Action_SymmetricRmsd']

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
