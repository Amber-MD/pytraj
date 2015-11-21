"""python write_action_pxd.py cpptraj_header_file.h"""
import re
from glob import glob

pxdlist = glob("Action_*.pxd")
print(pxdlist)

for action in pxdlist:
    # for action in [pxdlist[0],]:
    with open(action, 'r') as f0:
        txt = f0.read()

    oldline = "from Action cimport *"
    newline = "from pycpptraj.actions.Action cimport _Action, Action"
    txt = txt.replace(oldline, newline)

    fname = "bk/" + action
    with open(fname, 'w') as fh:
        fh.writelines(txt)
