# python3 ./enum_to_dict.py DataSet.h
# python3 ./enum_to_dict.py DataSet.h -se
# (-se = string to enum)
# (c) 2014 Hai Nguyen
from __future__ import print_function, absolute_import
import os
import CppHeaderParser
from util import print_blank_line, Line_codegen
from util import find_class
import sys

print(sys.argv[:])
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

for classname in list(cpp.classes.keys()):
    # create enum
    if cpp.classes[classname]['enums']['public']:
        for enumlist in cpp.classes[classname]['enums']['public']:
            print(indent + "# %s" % sys.argv[1])
            edictname = enumlist['name']
            if 'Type' in edictname:
                edictname = edictname.replace("Type", "")
            edictname = edictname + "Dict"
            print(edictname + " = {")

            # print key: value
            for enumvar in enumlist['values']:
                enumvarname = enumvar['name']
                if enum_to_string:
                    print(indent * 2 + '%s : "%s",' % (enumvarname, enumvarname))
                else:
                    print(indent * 2 + '"%s" : %s,' % (enumvarname, enumvarname))

            # end dict
            print("}")
