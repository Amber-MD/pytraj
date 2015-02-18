# (c) 2014 Hai Nguyen
import os
import CppHeaderParser
from util import print_blank_line, Line_codegen
from util import find_class
import sys

print sys.argv[:]
cpptrajsrc = os.environ["CPPTRAJHOME"] + "/src/"
filename = cpptrajsrc + sys.argv[1]
indent = " " * 4
classlist = find_class(cpptrajsrc)
cpp = CppHeaderParser.CppHeader(filename)

# format: PARM : "PARM",
enum_to_string = True
try:
    if sys.argv[2] and sys.argv[2] == '-se':
        # 'se' = 'string_to_enum'
        enum_to_string = False
except:
    pass

for classname in cpp.classes.keys():
    #create enum
    if cpp.classes[classname]['enums']['public']:
        for enumlist in cpp.classes[classname]['enums']['public']:
            print indent + "# %s" % sys.argv[1]
            enumname = enumlist['name']
            if 'Type' in edictname:
                edictname = edictname.replace("Type","")
            edictname = enumname + "Dict"
            print edictname + " = {"

            # print key: value
            for enumvar in enumlist['values']:
                enumvarname = enumvar['name']
                if enum_to_string:
                    print indent * 2 + '%s : "%s",' % (enumvarname, enumvarname)
                else:
                    print indent * 2 + '"%s" : %s,' % (enumvarname, enumvarname)

            # end dict
            print "}"

