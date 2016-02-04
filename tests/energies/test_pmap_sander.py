from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq

try:
    import sander
    has_sander = True
except ImportError:
    has_sander = False

code_global = '''
mm_options = sander.gas_input(8)
'''


@unittest.skipIf(not has_sander, 'skip if not having sander')
class TestSanderPmap(unittest.TestCase):

    def test_sander_pmap_simple(self):
        traj = pt.iterload('./data/Tc5b.x', './data/Tc5b.top')
        fname = traj.top.filename
        serial = pt.energy_decomposition(traj, prmtop=fname)['dihedral']
        parallel = pt.pmap(n_cores=4,
                           func=pt.energy_decomposition,
                           traj=traj,
                           prmtop=fname)['dihedral']
        aa_eq(serial, parallel)

    def test_sander_pmap_with_options(self):
        '''need to write mm_options as text
        '''
        code_local = '''
        mm_options = sander.gas_input(8)
        '''

        traj = pt.iterload('./data/Tc5b.x', './data/Tc5b.top')

        for code in [code_global, code_local]:
            data_parallel = pt.pmap(pt.energy_decomposition,
                                    traj,
                                    mm_options=code,
                                    n_cores=3,
                                    dtype='dict')

            mm_options = sander.gas_input(8)
            data_serial = pt.energy_decomposition(traj,
                                                  mm_options=mm_options,
                                                  dtype='dict')
            aa_eq(
                pt.tools.dict_to_ndarray(data_parallel),
                pt.tools.dict_to_ndarray(data_serial))


if __name__ == "__main__":
    unittest.main()
