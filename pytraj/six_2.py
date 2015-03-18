try:
    set = set
except NameError:
    from sets import Set as set
    set = set

try:
    from itertools import izip 
except ImportError:
    izip = zip
