from __future__ import print_function
from pytraj import io as mdio

@profile
def test_0():
    traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
    fa = traj[:]
    fa._fast_strip_atoms("@CA")
    fa._fast_strip_atoms("@CA")
    fa._fast_strip_atoms("@CA")

if __name__ == "__main__":
    test_0()
