import os
from time import time
from array import array
from pytraj import ArgList
from pytraj.cast_dataset import cast_dataset
from pytraj.actions.Action_Radgyr import Action_Radgyr
from pytraj.actions.Action_Molsurf import Action_Molsurf
from pytraj.actions.Action_Matrix import Action_Matrix
from pytraj.actions.Action_Strip import Action_Strip
from pytraj.actions.Action_Matrix import Action_Matrix
from pytraj.actions.Action_ClusterDihedral import Action_ClusterDihedral
from pytraj.CpptrajState import CpptrajState
from test_API.TestAPI import create_state, do_calculation
from pytraj.FrameArray import FrameArray

testdir = "../tests/data/"
topfile = testdir + "Tc5b.top"
trajinfile = testdir + "md1_prod.Tc5b.x"

state = create_state(top=topfile, trajin=trajinfile, ref=None)
#state.add_action(Action_Molsurf(), ArgList("surf"))
#state.add_action(Action_Molsurf(), ArgList("surf @CA"))
state.add_action(Action_Strip(), ArgList("strip !@CA"))
t0 = time()
state.run()
print("stripping time, ", time() - t0)
print(state.toplist[0].n_atoms)
#print state.datasetlist[0]
