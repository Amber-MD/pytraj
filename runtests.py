#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import time

# from: http://ascii.co.uk/art/batman
art = r'''
                |    |              _.-7
                |\.-.|             ( ,(_
                | a a|              \\  \,
                ) ["||          _.--' \  \\
             .-'  '-''-..____.-'    ___)  )\
            F   _/-``-.__;-.-.--`--' . .' \_L_
           |   l  {~~} ,_\  '.'.      ` __.' )\
           (    -.;___,;  | '- _       :__.'( /
           | -.__ _/_.'.-'      '-._ .'      \\
           |     .'   |  -- _                 '\,
           |  \ /--,--{ .    '---.__.       .'  .'
           J  ;/ __;__]. '.-.            .-' )_/
           J  (-.     '\'. '. '-._.-.-'--._ /
           |  |  '. .' | \'. '.    ._       \
           |   \   T   |  \  '. '._  '-._    '.
           F   J   |   |  '.    .  '._   '-,_.--`
           F   \   \   F .  \    '.   '.  /
          J     \  |  J   \  '.   '.    '/
          J      '.L__|    .   \    '    |
          |   .    \  |     \   '.   '. /
          |    '    '.|      |    ,-.  (
          F   | ' ___  ',._   .  /   '. \
          F   (.'`|| (-._\ '.  \-      '-\
          \ .-'  ( L `._ '\ '._ (
     snd  /'  |  /  '-._\      ''\
              `-'
'''

my_script = sys.argv[0]

try:
    need_help = sys.argv[1] in ['help', '-help', '--help']
except IndexError:
    need_help = False 
try:
    do_simple_test = sys.argv[1] in ['simple', 'minimal', '-simple',
                                     '-minimal', 'sim']
except:
    do_simple_test = False

if need_help:
    print("Usage:")
    print("    short testing: python %s simple" % my_script)
    print("    long testing: python %s" % my_script)
    print("Note: long testing requires nose and coverage, which are easily installed by `pip install`")
    sys.exit(0)

print("start testing. Go to ./tests folder")
os.chdir("./tests/")

if do_simple_test:
    os.system("python ./run_simple_test.py")
    print('\nHAPPY COMPUTING')
    print(art)
    sys.exit(0)
else:
    os.system("nosetests --with-coverage --cover-package pytraj -vs .")

print('\nHAPPY COMPUTING')
print(art)
