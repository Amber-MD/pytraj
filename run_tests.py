#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import subprocess
import argparse

bin = sys.prefix + '/bin/'
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
parser = argparse.ArgumentParser(
    description='run test. Full test requires nose and coverage packages')
parser.add_argument('-s', '--simple', action='store_true', help='quick run')
parser.add_argument(
    '-c', '--with-coverage', action='store_true', help='coverage report')
parser.add_argument('-n', '--ncores', default=1, help='n_cores')
args = parser.parse_args()

print(bin)
print("start testing. Go to ./tests folder")
os.chdir("./tests/")

if args.simple:
    print('running minimal test\n')
    subprocess.check_call(
        "{bin}/python ./run_simple_test.py".format(bin=bin), shell=True)
else:
    print('running full test\n')
    if args.with_coverage:
        print('with coverage')
        cm_list = "{bin}/pytest --cov=pytraj --cov-report=html -vs -rsx .".format(
            bin=bin).split()
    else:
        print('without coverage\n')
        cm_list = "{bin}/pytest -vs .".format(bin=bin).split()
    if int(args.ncores) > 1:
        cm_list.extend(['-n', args.ncores])
    print("using {} cores".format(args.ncores))
    subprocess.check_call(cm_list)

print('\nHAPPY COMPUTING')
print(art)
