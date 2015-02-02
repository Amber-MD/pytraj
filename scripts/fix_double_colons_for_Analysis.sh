#!/bin/sh

python ../scripts_for_wrapping/codegen_pxd.py $1.h | sed "s/_Analysis::RetType/RetType/g"
