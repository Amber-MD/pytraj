#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles nastruct.in BP.*.dat BPstep.*.dat bases.pdb baseaxes.pdb basepairaxes.pdb Helix.*.dat Param.pdb

# Test 2
INPUT="-i nastruct.in"
cat > nastruct.in <<EOF
parm ../adh026.3.pdb
trajin ../adh026.3.pdb 
nastruct naout adh026.dat
EOF
RunCpptraj "NAstruct command test."
DoTest BP.adh026.dat.save BP.adh026.dat
DoTest BPstep.adh026.dat.save BPstep.adh026.dat
DoTest Helix.adh026.dat.save Helix.adh026.dat
CheckTest

EndTest

exit 0
