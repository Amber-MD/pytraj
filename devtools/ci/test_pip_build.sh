#!/bin/sh

./devtools/mkrelease
(cd dist && python ../scripts/build_wheel.py *gz)
