#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import time

bin = sys.prefix + '/bin/'
print(bin)

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
    sys.argv.remove('--with-coverage')
    with_coverage = True
except ValueError:
    with_coverage = False

try:
    sys.argv.remove('x')
    with_coverage = True
except ValueError:
    pass

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
    print("    short testing: python {} simple".format(my_script))
    print("    long testing: python {}".format(my_script))
    print(
        "    long testing with code coverage: python {} --with-coverage".format(
            my_script))
    print(
        "Note: long testing requires nose and coverage, which are easily installed by `pip install`")
    sys.exit(1)

print("start testing. Go to ./tests folder")
os.chdir("./tests/")

if do_simple_test:
    os.system("python ./run_simple_test.py")
    print('\nHAPPY COMPUTING')
    print(art)
    sys.exit(1)
else:
    if with_coverage:
        os.system("{bin}/nosetests --with-coverage --cover-package pytraj -vs .".format(bin=bin))
    else:
        os.system("{bin}/nosetests -vs .".format(bin=bin))

print('\nHAPPY COMPUTING')
print(art)
