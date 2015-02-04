from __future__ import print_function
from pytraj.actions import allactions

adict = {}

for key in allactions.__dict__.keys():
    if "Action_" in key:
        act = key.split('Action_')[1]
        #print (act)
        # create dict of action objects
        adict[act.lower()] = allactions.__dict__[key]

# make another dict to convert something like 'MolSurf' to 'molsurf'
