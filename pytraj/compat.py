from __future__ import absolute_import

from .externals.six.moves import range, map
from .externals.six import string_types, callable

try:
    set = set
except NameError:
    from sets import Set as set
    set = set

try:
    from itertools import izip 
except ImportError:
    izip = zip

zip = izip
