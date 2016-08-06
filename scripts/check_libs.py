#!/usr/bin/env python
#

import warnings

warnings.filterwarnings('ignore', category=UserWarning)

suggested_modules = ['sphinx', 'sander', 'parmed', 'numpy', 'mpi4py',
                     'sklearn', 'pandas', 'matplotlib', 'seaborn',
                     'sphinx_bootstrap_theme']

suggested_modules.remove('sander')

for module in suggested_modules:
    try:
        lib = __import__(module)
    except ImportError as e:
        print('\ncheck README.md\n')
        print('or use scripts/install_all.sh to install all required packages')
        raise ImportError('require {}'.format(module))
