#!/bin/sh

#git subtree pull --prefix AmberTools/src/pytraj pytraj_github master

# merge all commmits to single one
git pull --squash -s recursive -X subtree=AmberTools/src/pytraj pytraj_github master

# need to git commit after squashing

# log
# 2030  git pull --squash -s recursive -X subtree=AmberTools/src/pytraj pytraj_github master
# 2031  git status
# 2032  git status | less
# 2033  git commit
# 2034  git stauts
# 2035  git status
# 2036  git status -uno
# 2037  git log 
# 2038  git checkout 94e628c0c862d86ac720760dc23dd537ffab1c46
# 2039  ls
# 2040  ls AmberTools/
# 2041  ls AmberTools/src/
# 2042  ls AmberTools/src/sander/
# 2043  ls AmberTools/src/sander/*.F90
# 2044  ls
# 2045  git checkout master
