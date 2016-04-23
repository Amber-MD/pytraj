#!/usr/bin/env python
#

import warnings

warnings.filterwarnings('ignore', category=UserWarning)

for module in ['sphinx', 'sander', 'parmed', 'numpy', 'mpi4py',
               'pandas', 'matplotlib', 'seaborn', 'sphinx_bootstrap_theme']:
    try:
        lib = __import__(module)
    except ImportError as e:
        print('\ncheck README.md\n')
        print('or use scripts/install_all.sh to install all required packages')
        raise ImportError('require {}'.format(module))
