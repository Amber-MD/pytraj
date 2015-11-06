#!/bin/sh

# adpated from ParmEd
# https://github.com/ParmEd/ParmEd

# This script runs pyflakes and filters out all of the known violations to
# provide an output of what the linters find.

which pyflakes > /dev/null

if [ $? -ne 0 ]; then
    echo "pyflakes not found... cannot run check."
    exit 1
fi

notfound() {
    echo "pytraj package directory not found!"
    exit 1
}

test -d pytraj || notfound

pyflakes pytraj 2>&1 | \
    grep -v -e "__init__.py" \
            -e "pytraj/externals/six.py" \
            -e "pytraj/run_tests.py" \
            -e "^parmed\/unit" \
            -e "^parmed\/utils\/six.py" \
            -e "^parmed\/utils\/fortranformat" | tee pyflakes.log

nfail=`cat pyflakes.log | wc -l`

if [ $nfail -gt 0 ]; then
    echo "Detected pyflakes failures"
    exit 1
fi

echo "pyflakes reported clean exit"
/bin/rm -f pyflakes.log
exit 0
