#!/bin/sh

python setup.py register
python setup.py sdist upload -r pypi
