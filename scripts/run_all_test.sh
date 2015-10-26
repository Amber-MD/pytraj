#!/bin/sh

for test in test*.py; do
echo "File: " $test
python $test
done
