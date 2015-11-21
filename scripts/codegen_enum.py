# (c) 2014 - Hai Nguyen
import os
import sys
import CppHeaderParser
from util import find_class

'''python ./codegen_enum.py dict
'''

CPPTRAJSRC = os.environ['CPPTRAJHOME'] + "/src/"


def create_enum_of_dict(fname, mode='', cpptrajsrc=CPPTRAJSRC):
    #cpptrajsrc = os.environ['AMBERHOME'] + "AmberTools/src/cpptraj/src/"
    fname_full = cpptrajsrc + fname

    indent = " " * 4
    tmpindent = " " * 4
    classlist = find_class(cpptrajsrc)
    cpp = CppHeaderParser.CppHeader(fname_full)

    make_dict = False
    if mode == 'dict':
        make_dict = True
        tmpindent = ""
    else:
        make_dict = False

    # make assumption that there's only one class in header file
    #classname = list(cpp.classes.keys())[0]
    # use -1 for MetaData, else 0
    classname = list(cpp.classes.keys())[-1]
    print(cpp.classes[classname]['enums'])
    if cpp.classes[classname]['enums']['public']:
        print(classname)
        for enumlist in cpp.classes[classname]['enums']['public']:
            print("\n")
            if make_dict:
                print("from %s cimport *" % fname)
            else:
                print('cdef extern from "%s":' % fname)
            enumname = enumlist['name']
            enumext = classname + "::" + enumname
            if not make_dict:
                print(indent + 'ctypedef enum %s "%s":' % (enumname, enumext))
            else:
                enumname = enumname.replace("Type", "")
                print("%sDict = {" % enumname)
            for enumvar in enumlist['values']:
                enumvarname = enumvar['name']
                enumvarnameext = classname + "::" + enumvarname
                if not make_dict:
                    print(indent * 2 + '%s "%s"' % (enumvarname, enumvarnameext))
                else:
                    print(indent + '"%s" : %s, ' % (enumvarname, enumvarname))
            if make_dict:
                print(indent + "}")

if __name__ == '__main__':
    fname = sys.argv[1]
    try:
        mode = sys.argv[2]
    except IndexError:
        mode = ""
    create_enum_of_dict(fname, mode=mode)
