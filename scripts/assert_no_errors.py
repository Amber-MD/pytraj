#!/usr/bin/env python

from __future__ import print_function

# try to find html pages that having compiling error
# (such as ImportError)

from glob import glob

excluded_files = set(['latest/faq.html',])
error_files = []

for fn in glob('latest/*html') + glob('latest/*/*html'):
    with open(fn) as fh:
        content = fh.read()
        if 'ImportError' in content:
            error_files.append(fn)
        if 'DeprecationWarning:' in content:
            error_files.append(fn)
        if 'NameError:' in content:
            error_files.append(fn)
        if 'OSError' in content:
            error_files.append(fn)
        if 'NotImplementedError' in content:
            error_files.append(fn)
        if 'AttributeError' in content:
            error_files.append(fn)

print(set(error_files))
assert set(error_files) - excluded_files == set()
