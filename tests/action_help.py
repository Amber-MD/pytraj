#!/usr/bin/env python
import sys
from pytraj import adict

key = sys.argv[1].lower()
adict[key].help()
