# (c) 2014 Hai Nguyen

# TODO: should re-write, looks ugly
import os
import CppHeaderParser
from util import print_blank_line, Line_codegen
from util import find_class
import sys

#cpptrajsrc = os.environ['AMBERHOME'] + "/AmberTools/src/cpptraj/src/"
cpptrajsrc = os.environ['CPPTRAJHOME'] + "/src/"
filename = cpptrajsrc + sys.argv[1]
short_filename = filename.split("/")[-1]
indent = " " * 4
classlist = find_class(cpptrajsrc)
cpp = CppHeaderParser.CppHeader(filename)

# check if need extract line for Action_*.pyx or Analysis_*.pyx classes
need_extra_line = False
if short_filename.startswith("Action") or short_filename.startswith("Analysis"):
    need_extra_line = True
    # get Action's name from short_filename
    # Action_Rmsd.h --> Action_Rmsd
    action = short_filename.split(".")[0]
    # actionroot = "Action" or "Analysis"?
    actionroot = action.split("_")[0]

# print header line "c++" so Cython know it is c++ code
# (adding to setup.py seems not work)
print("# distutils: language = c++")

# add vector, string if having ones
stdlist = ['vector', 'string']
with open(filename, 'r') as fh:
    linelist = " ".join(fh.readlines())
    for word in stdlist:
        if word in linelist:
            print("from libcpp.%s cimport %s" % (word, word))

for f_include in cpp.includes:
    # remove "
    if f_include.startswith('"'):
        f_include = f_include.split('"')[1]
    # import stuff (in cpptraj header files, need to add libcpp.* too)
    if not f_include.startswith("<"):
        print("from %s cimport *" % (f_include.split(".")[0]))

print_blank_line(2)
print('cdef extern from "%s": ' % short_filename)

# make assumption that there's only one class in header file
for classname in list(cpp.classes.keys()):
    # create enum
    if cpp.classes[classname]['enums']['public']:
        for enumlist in cpp.classes[classname]['enums']['public']:
            print(indent + "# %s" % sys.argv[1])
            enumname = enumlist['name']
            enumext = classname + "::" + enumname
            print(indent + 'ctypedef enum %s "%s":' % (enumname, enumext))
            for enumvar in enumlist['values']:
                enumvarname = enumvar['name']
                enumvarnameext = classname + "::" + enumvarname
                print(indent * 2 + '%s "%s"' % (enumvarname, enumvarnameext))

    # declare cpp class
    extcl = "_" + classname
    if not need_extra_line:
        print('%scdef cppclass %s "%s":' % (indent, extcl, classname))
    else:
        print('%scdef cppclass %s "%s" (_%s):' % (indent, extcl, classname, actionroot))
    methods = cpp.classes[classname]['methods']['public']
    for method in methods:
        line = Line_codegen(method['debug'])
        # move to line's methods?
        if 'virtual' in line.myline:
            line.add_sharp()
        line.remove_std_namespace()
        line.swap_const()
        line.add_under_score_to_class(classlist)
        line.replace_others()
        line.swap_const()
        line.remove_unsupported()
        line.remove_preassignment()
        print(indent * 2 + line.myline)
    print_blank_line(2)

for classname in list(cpp.classes.keys()):
    if not need_extra_line:
        print("cdef class %s:" % classname)
    else:
        print("cdef class %s (%s):" % (classname, actionroot))
    print("%scdef _%s* thisptr" % (indent, classname))
    print()
