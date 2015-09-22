#!/usr/bin/env python
"""Get help: python ./pytraj.py --help"""

import sys
from pytraj.cpptraj import Cpptraj

print(sys.argv[:])
Cpptraj().run([x.encode() for x in sys.argv[:]])
