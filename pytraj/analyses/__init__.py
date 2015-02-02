""""""
from __future__ import absolute_import
from glob import glob
from pytraj.analyses.allanalyses import *
from pytraj.analyses import allanalyses

# get analysislist from `allanalyses module`
analysislist = []
for act in allanalyses.__dict__.keys():
    if 'Analysis' in act:
        analysislist.append(act)

__all__ = analysislist + ['allanalyses', ]
