#!/usr/bin/env python

'''
code generation for Action_*.pyx

CPPTRAJHOME must be set
'''
import sys
import os

cpptrajhome = sys.argv[1]

if not cpptrajhome:
    raise EnvironmentError('must set CPPTRAJHOME')

actionlist = []
cpptraj_analist = []

action_pyx = os.path.join(
        os.path.dirname(__file__), '../..',
        'pytraj/analysis/c_action/c_action.pyx')


def get_header(action_pyx):
    with open(action_pyx) as fh:
        lines = fh.readlines()
        for index, line in enumerate(lines):
            if line.startswith('cdef class Action_Angle'):
                break
        return ''.join(lines[:index])


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

text_pxd = """
cdef extern from "{action_name}.h": 
    cdef cppclass _{action_name} "{action_name}" (_Action) nogil:
        _{action_name}()
        _DispatchObject * Alloc() 
        void Help() 


cdef class {action_name}(Action):
    cdef _{action_name}* thisptr
"""

with open('tmp.pyx', 'w') as fh:
    fh.write(get_header(action_pyx))
    for action in actionlist:
        tmp = text.format(action_name=action)
        # tmp = text_pxd.format(action_name=action)
        fh.write(tmp)

print(action_pyx)
