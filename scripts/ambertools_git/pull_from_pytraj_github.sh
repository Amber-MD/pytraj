#!/bin/sh

git fetch pytraj_github
# merge all commmits to single one
git pull -s recursive -X subtree=AmberTools/src/pytraj -X theirs --squash pytraj_github master
