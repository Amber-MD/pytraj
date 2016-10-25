#!/usr/bin/env python
#

import warnings

warnings.filterwarnings('ignore', category=UserWarning)

amber_modules = ['parmed']
suggested_modules = ['sphinx', 'numpy', 'mpi4py',
                     'sklearn', 'pandas', 'matplotlib', 'seaborn',
                     'sphinx_bootstrap_theme',
                     'IPython',
                     'notebook',
                     'runipy']
suggested_modules += amber_modules

missing_packages = []
for module in suggested_modules:
    try:
        lib = __import__(module)
    except ImportError as e:
        missing_packages.append(module)

command = '    conda install {}'.format(' '.join(missing_packages))

if 'sphinx_bootstrap_theme' in command:
    command = command.replace('sphinx_bootstrap_theme', '')

if missing_packages:
    print('missing: ')
    for module in missing_packages:
        print(module)

    print('\nTRY: ')
    print(command)
    if 'sphinx_bootstrap_theme' in missing_packages:
        print('pip install sphinx_bootstrap_theme')
    print("")
    raise ImportError()
