#!/bin/sh

fname=$1

./grep_wrap.sh "$fname::" $fname.cpp | sed "/\/\//d"
