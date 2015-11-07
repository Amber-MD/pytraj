""""""
from __future__ import absolute_import
from pytraj.actions import CpptrajActions

actionlist = []
for act in CpptrajActions.__dict__.keys():
    if 'Action' in act:
        actionlist.append(act)

__all__ = actionlist
__doc__ = "\n".join(__all__)
