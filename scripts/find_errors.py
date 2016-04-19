#!/usr/bin/env python

# try to find html pages that having compiling error
# (such as ImportError)

from glob import glob

for fn in glob('latest/*html') + glob('latest/*/*html'):
    if 'ImportError' in open(fn).read():
        print(fn)
