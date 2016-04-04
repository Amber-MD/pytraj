#!/usr/bin/env python

'''
code generation for Action_*.pyx

CPPTRAJHOME must be set
'''
import os

cpptrajhome = os.environ.get('CPPTRAJHOME', '')

if not cpptrajhome:
    raise EnvironmentError('must set CPPTRAJHOME')

actionlist = []
cpptraj_analist = []

with open(cpptrajhome + '/src/Command.cpp') as fh:
    lines = fh.readlines()
    actionlist += [line.split('\n')[0].split()[1].replace('"', '').replace(',', '').replace('.h', '')
                        for line in lines if line.startswith('#include "Action_')]
exlucdedList = ['Action_CreateReservoir',]

for excluded_action in exlucdedList:
    actionlist.remove(excluded_action)

text = """
cdef class {action_name}(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _{action_name}()
        self.thisptr = <_{action_name}*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()
"""

for action in actionlist:
    tmp = text.format(action_name=action)
    print(tmp)
