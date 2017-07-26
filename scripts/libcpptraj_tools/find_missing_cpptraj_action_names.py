#!/usr/bin/env python
'''
Example
-------
$ python ./find_missing_cpptraj_action_names.py
{'Action_CreateReservoir', 'Action_Unstrip'}
'''

import os
import sys

cpptrajhome = sys.argv[1]

cpptraj_actlist = []
cpptraj_analist = []

with open(cpptrajhome + '/src/Command.cpp') as fh:
    lines = fh.readlines()
    cpptraj_actlist += [line.split('\n')[0].split()[1].replace('"', '').replace('.h', '') for line in lines if '#include "Action_' in line]
    cpptraj_analist += [line.split('\n')[0].split()[1].replace('"', '').replace('.h', '') for line in lines if '#include "Analysis_' in line]

from pytraj.analysis.c_action import c_action as CA
from pytraj.analysis.c_analysis import c_analysis as CAnal
pytraj_actlist = [key for key in dir(CA) if key.startswith("Action_")]
pytraj_analist = [key for key in dir(CAnal) if key.startswith("Analysis_")]

print(set(pytraj_actlist) ^ set(cpptraj_actlist))
print(set(pytraj_analist) ^ set(cpptraj_analist))
