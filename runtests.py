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
    verbose = sys.argv[1] in ['-verbose', 'verbose', '-v']
except:
    verbose = False
try:
    need_help = sys.argv[1] in ['help', '-help', '--help']
except:
    need_help = False 
try:
    do_simple_test = sys.argv[1] in ['simple', 'minimal', '-simple',
                                     '-minimal', 'sim']
except:
    do_simple_test = False

if need_help:
    print("Usage:")
    print("    short testing: python %s simple" % my_script)
    print("    verbose run (long): python %s verbose" % my_script)
    print("    quiet run (long): python %s" % my_script)
    sys.exit(0)

print("start testing. Go to ./tests folder")
os.chdir("./tests/")

if do_simple_test:
    os.system("python ./run_simple_test.py")
    print('\nHAPPY COMPUTING')
    print(art)
    sys.exit(0)
if verbose:
    os.system('python get_unittest_files.py')
    os.system('sh TestListTravis.sh')
else:
    print('quite run')
    os.system("python ./run_all_and_find_fails.py")

print('\nHAPPY COMPUTING')
print(art)
