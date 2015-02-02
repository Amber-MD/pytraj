""""""
from __future__ import absolute_import
from glob import glob
from pytraj.actions.allactions import *
from pytraj.actions import allactions

# get actionlist from `allactions module`
#actionlist = [act if ('Action' in act) for act in allactions.__dict__.keys()]
actionlist = []
for act in allactions.__dict__.keys():
    if 'Action' in act:
        actionlist.append(act)

__all__ = actionlist + ['allactions',]
__doc__ = "\n".join(__all__)
