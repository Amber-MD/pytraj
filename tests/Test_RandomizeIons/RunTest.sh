#!/bin/bash

. ../MasterTest.sh

CleanFiles ptraj.in random.crd 
CheckZlib
INPUT="ptraj.in"
TOP="adh206.ff10.tip3p.parm7.gz"
cat > ptraj.in <<EOF
trajin adh206.tip3p.rst7.gz
randomizeions @Na+ around :1-16 by 5.0 overlap 3.0
trajout random.crd title "Test"
EOF
RunCpptraj "randomizeions test"
DoTest random.crd.save random.crd
CheckTest

EndTest

exit 0
