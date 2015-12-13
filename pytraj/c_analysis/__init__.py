""""""
from __future__ import absolute_import
from glob import glob
from pytraj.analyses import CpptrajAnalyses

analysislist = []
for act in CpptrajAnalyses.__dict__.keys():
    if 'Analysis' in act:
        analysislist.append(act)

__all__ = analysislist
