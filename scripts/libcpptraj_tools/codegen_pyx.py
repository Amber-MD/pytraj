# (c) 2014 Hai Nguyen
# TODO: should re-write, looks ugly
# Note for me: method's name (example)
# cppframe['methods']['public'][12]['name']
import os
import CppHeaderParser
from util import print_blank_line, Line_codegen
from util import find_class
import sys

cpptrajsrc = os.environ['AMBERHOME'] + "/AmberTools/src/cpptraj/src/"
#cpptrajsrc = os.environ['CPPTRAJHOME'] + "/src/"
file = cpptrajsrc + sys.argv[1]
# try:
#    # virtual or not?
#    option = sys.argv[2]
#    if option == 'virtual':
#        include_virtual = True
# except:
#    include_virtual = False
include_virtual = True

indent = " " * 4
classlist = find_class(cpptrajsrc)
cpp = CppHeaderParser.CppHeader(file)

# print header line "c++" so Cython know it is c++ code
# (adding to setup.py seems not work)
print("# distutils: language = c++")
print_blank_line(2)

# make assumption that there's only one class in header file
for classname in list(cpp.classes.keys()):
    # declare cpp class
    print('cdef class %s:' % (classname))
    print('    def __cinit__(self):')
    print('        self.thisptr = new _%s()' % classname)
    print_blank_line(1)
    print('    def __dealloc__(self):')
    print('        del self.thisptr')
    print_blank_line(1)

    methods = cpp.classes[classname]['methods']['public']
    for method in methods:
        line = Line_codegen(method['debug'])
        # if 'virtual' not in line.myline:
        if include_virtual:
            line.remove_std_namespace()
            # line.remove_unsupported()
            line.swap_const()
            line.replace_others()
            # call swap_const() again to change "vector[int] const& to const vector[int]&"
            line.swap_const()
            line.remove_word()
            line.remove_unsupported()
            line.remove_preassignment()
            line.insert_self_word()
            if not line.has_ignored_words():
                print("%sdef %s" % (indent, line.myline))
                print_blank_line(1)
