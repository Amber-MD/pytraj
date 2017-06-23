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
        possible_errors = [
                'ImportError',
                'TypeError',
                'DeprecationWarning',
                'NameError',
                'OSError',
                'NotImplementedError',
                'AttributeError'
        ]
        for error in possible_errors:
            if error in content:
                error_files.append(fn)

print(set(error_files))
assert set(error_files) - excluded_files == set()
