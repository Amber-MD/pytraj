import os
from time import time
from array import array
from pytraj import ArgList
from pytraj.datasets import cast_dataset
from pytraj.actions.CpptrajActions import Action_Radgyr
from pytraj.actions.CpptrajActions import Action_Molsurf
from pytraj.actions.CpptrajActions import Action_Matrix
from pytraj.actions.CpptrajActions import Action_Strip
from pytraj.actions.CpptrajActions import Action_Matrix
from pytraj.actions.CpptrajActions import Action_ClusterDihedral
from pytraj.core.CpptrajState import CpptrajState
from test_API.TestAPI import create_state, do_calculation
from pytraj.Trajectory import Trajectory

#testdir = os.environ['PYCPPTRAJ_HOME'] + "/tests/Cpptraj_test/"
#topfile = testdir + "Test_Matrix/1rrb_vac.prmtop"
#trajinfile = testdir + "Test_Matrix/1rrb_vac.mdcrd"
testdir = "./data/"
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
