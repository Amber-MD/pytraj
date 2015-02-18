import os
import pytraj
from importlib import import_module

# get action list
pylist = "../pyxlist.txt"

actions = []
with open(pylist, 'r') as fh:
    lines = fh.readlines()
    for line in lines:
        if line.startswith("actions/"):
            try:
                actions.append(line.split("/")[-1].split()[0])
            except: pass

removedList = ["Action", "ActionFrameCounter", "Action_CreateReservoir"]
for removedAction in removedList:
    actions.remove(removedAction)

# get Help from actionlist
def get_help():
    for act_name in actions:
        module = '.actions.' + act_name
        if "Action_" in module:
            mod = import_module(module, package='pytraj')
            classname = mod.__getattribute__(act_name)
            print("Help for %s" % act_name)
            classname().help()
            print()

get_help()
