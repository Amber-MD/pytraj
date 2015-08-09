import numpy as np
from pytraj.actions.CpptrajActions import Action_Distance
from test_API.TestAPI import create_state, do_calculation

state = create_state(
    top="./data/Tc5b.top",
    trajin="./data/md1_prod.Tc5b.x",
    ref=None)
distance = do_calculation(
    action=Action_Distance(),
    command="distance :2@CA :10@CA",
    state=state)

distance[:10]

# make sure to reproduce cpptraj output
cppout = np.loadtxt(
    "./data/CAres2_CAres10.Tc5b.dat",
    skiprows=1).transpose()[1]
print(cppout[:10])

dsize = len(distance)
np.testing.assert_almost_equal(distance, cppout[:dsize], decimal=3)
print("Kool")
