from __future__ import print_function
from pytraj.testing import test_if_having, aa_eq
from pytraj import io
import pytraj.common_actions as pyca

# See Also: ../tests/test_calc_energies.py


@test_if_having("sander")
@test_if_having("parmed")
def main():
    # require 'sanderapi' and 'parmed' packages in AmberTools15
    # http://ambermd.org/#AmberTools
    # more info about `sanderapi` in Amber15 manual
    # http://ambermd.org/doc12/Amber15.pdf (page 341)
    traj = io.iterload("../tests/data/Tc5b.x",
                       "../tests/data/Tc5b.top")
    print(traj)
    edict0 = pyca.energy_decomposition(parm="../tests/data/Tc5b.top",
                                       traj=traj,
                                       igb=8,
                                       input_options=None)
    print(edict0)

    # edict1: use default igb=8
    edict1 = pyca.energy_decomposition(traj=traj)
    print(edict1)

    aa_eq(edict0['tot'], edict1['tot'])

    # update `input` and get minimal output
    import sander
    inp = sander.gas_input(6)
    edict2 = pyca.energy_decomposition(traj, input_options=inp, mode='minimal')
    print(edict2)


if __name__ == "__main__":
    main()
