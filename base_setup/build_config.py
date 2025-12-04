"""
Build configuration management - centralized flag and option handling
"""
import os
import sys
import builtins

# Import from local base_setup module
from .base_setup import get_version_info, CleanCommand

class BuildConfig:
    """Centralized configuration for build process"""

    def __init__(self):
        # Parse command line flags
        self.amber_release = self._check_flag('--amber_release')
        self.disable_openmp = self._check_flag('--disable-openmp')
        self.use_amberlib = not self._check_flag('--disable-amberlib')
        self.use_prebuilt = self._check_flag('--use-pre-cythonized')
        self.debug = self._check_flag('-debug')
        self.cythonize_only = self._check_flag('--cythonize')

        # Determine compilation mode
        self.compile_extensions = self._should_compile()
        self.is_sdist = 'sdist' in sys.argv
        self.use_pip = self._detect_pip_install()

        # OpenMP flag
        self.openmp_flag = '-openmp' if not self.disable_openmp else ''

        # Version info
        self.version, self.git_revision = get_version_info()

        # Package configuration
        self.packages = self._get_packages()
        self.classifiers = self._get_classifiers()
        self.cmdclass = {'clean': CleanCommand}

        # Set the global flag for pytraj detection
        builtins.__PYTRAJ_SETUP__ = True

    def _check_flag(self, flag):
        """Check and remove flag from sys.argv"""
        try:
            sys.argv.remove(flag)
            return True
        except ValueError:
            return False

    def _should_compile(self):
        """Determine if C extensions should be compiled"""
        if '--help' in sys.argv or '-h' in sys.argv or '--help-commands' in sys.argv:
            return False
        return (self.use_pip or
                'install' in sys.argv or
                'build' in sys.argv or
                'build_ext' in sys.argv)

    def _detect_pip_install(self):
        """Check if being installed via pip"""
        return any(arg in sys.argv for arg in ['egg_info', 'pip', '--no-deps'])

    def _get_packages(self):
        """Return package list"""
        return [
            'pytraj', 'pytraj.utils', 'pytraj.builder',
            'pytraj.actions', 'pytraj.analysis',
            'pytraj.analysis.c_action', 'pytraj.analysis.c_analysis',
            'pytraj.datasets', 'pytraj.externals',
            'pytraj.trajectory', 'pytraj.trajectory.c_traj',
            'pytraj.topology', 'pytraj.datafiles',
            'pytraj.datafiles.ala3', 'pytraj.datafiles.tz2',
            'pytraj.datafiles.dpdp', 'pytraj.datafiles.trpcage',
            'pytraj.datafiles.remd_ala2', 'pytraj.math',
            'pytraj.core', 'pytraj.parallel', 'pytraj.cluster',
            'pytraj.visualization', 'pytraj.serialize',
            'pytraj.sandbox', 'pytraj.testing',
        ]

    def _get_classifiers(self):
        """Return classifier list"""
        return [
            'Development Status :: 5 - Production/Stable',
            'Operating System :: Unix',
            'Operating System :: MacOS',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: Python :: 3.10',
            'Programming Language :: Python :: 3.11',
            'Programming Language :: Python :: 3.12',
            'Programming Language :: Cython',
            'Programming Language :: C',
            'Programming Language :: C++',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry',
        ]

    def get_cython_directives(self):
        """Get Cython compilation directives"""
        base_directives = {
            'embedsignature': True,
            'boundscheck': False,
            'wraparound': False,
            'auto_pickle': False,
            'language_level': 3
        }

        if self.debug:
            base_directives.update({
                'profile': True,
                'linetrace': True,
                'binding': True
            })

        return base_directives

    def get_define_macros(self):
        """Get C preprocessor macros"""
        if self.debug:
            return [('CYTHON_TRACE', 1)]
        return []