#!/bin/sh

pytrajamber=$1
pytrajgit=$2
git diff --stat $pytrajamber:AmberTools/src/pytraj pytraj_github/$pytrajgit
