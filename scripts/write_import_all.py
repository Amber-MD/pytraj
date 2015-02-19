import re
from glob import glob

pxdlist = glob("Action*.pxd")
actionlist = []

for pxd in pxdlist:
    actionlist.append(pxd.split(".")[0])

exlucdedList = []
for excluded_action in exlucdedList:
    actionlist.remove(excluded_action)

for action in actionlist:
    print("from pycpptraj.actions.%s import %s" % (action, action))
