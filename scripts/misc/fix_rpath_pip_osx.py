#!/usr/bin/env python

from conda_build import post
from glob import glob

for fn in (glob('pytraj/*so') + glob('pytraj/*/*so')):
    post.mk_relative_osx(fn, 'pytraj')

import pytraj as pt
