import re
from glob import glob

pxdlist = glob("Action*.pxd")
actionlist = []

for pxd in pxdlist:
    try:
        action = re.findall("Action_(.+?).pxd", pxd)[0]
    except:
        pass
    actionlist.append(action)

exlucdedList = ['Rmsd', 'Dihedral']
for excluded_action in exlucdedList:
    actionlist.remove(excluded_action)

#print actionlist
text = """# distutils: language = c++
from cython.operator cimport dereference as deref


cdef class Action_ACTION_NAME (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_ACTION_NAME()
        self.thisptr = <_Action_ACTION_NAME*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        \"""return a function-pointer object to be used with ActionList class
        \"""
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()
""" 

for action in actionlist:
    tmp = text.replace("ACTION_NAME", action)
    fname = "Action_" + action + ".pyx"
    with open(fname, 'w') as fh:
        fh.writelines(tmp)
