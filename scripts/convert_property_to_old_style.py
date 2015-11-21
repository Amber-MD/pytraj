"""
convert:

@property
def x(self):
    return get_value((x)

to:
property x:
    def __get__(self):
        return get_value(x)

Why? Cython has not yet supported the new property decorator
"""
import sys
import re

indent = ' ' * 4
lines = open(sys.argv[1], 'r').readlines()
for i, line in enumerate(lines):
    if "@property" in line:
        block = lines[i + 1]
        fname = re.findall(r'def (.+?)\(self\)', block)[0]
        print(indent + 'property ' + fname + ':')
        print(indent * 2 + 'def __get__(self):')
        print(indent * 3 + 'return self.thisptr.' + fname + '()')
        print()
