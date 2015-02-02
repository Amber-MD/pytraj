#!/bin/sh

python setup.py register -r pypitest
python setup.py sdist upload -r pypitest
