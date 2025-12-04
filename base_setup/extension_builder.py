"""
C Extension building logic - handles Cython compilation and cpptraj integration
"""
import os
import sys
import subprocess
from glob import glob
from setuptools import Extension

# Import from local base_setup module
from .base_setup import (
    get_cpptraj_info, get_ext_modules, get_package_data,
    check_cython, try_updating_libcpptraj, add_openmp_flag,
    setenv_cc_cxx, write_version_py
)class ExtensionBuilder:
    """Handles building C extensions with cpptraj integration"""

    def __init__(self, config, platform):
        self.config = config
        self.platform = platform
        self.rootname = os.getcwd()
        self.pytraj_src = os.path.join(self.rootname, "pytraj")
        self.cpptraj_home = os.environ.get('CPPTRAJHOME', '')
        self.cpptraj_included = os.path.exists("./cpptraj/")
        self.pytraj_home = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

        # Write version file
        write_version_py()

        # Get cpptraj information
        self.cpptraj_info = self._get_cpptraj_info()

        # Setup Cython
        self.need_cython, self.cmdclass, self.cythonize = self._setup_cython()

        # Get library files
        self.libcpptraj_files = self._get_libcpptraj_files()

    def _get_cpptraj_info(self):
        """Get cpptraj installation information"""
        return get_cpptraj_info(
            rootname=self.rootname,
            cpptraj_home=self.cpptraj_home,
            cpptraj_included=self.cpptraj_included,
            compile_c_extension=self.config.compile_extensions,
            pytraj_home=self.pytraj_home,
            openmp_flag=self.config.openmp_flag,
            use_amberlib=self.config.use_amberlib
        )

    def _setup_cython(self):
        """Setup Cython compilation"""
        return check_cython(
            is_released=False,  # Based on original setup.py logic
            cmdclass=self.config.cmdclass,
            min_version='0.21',
            use_prebuilt_cythonized_files=self.config.use_prebuilt
        )

    def _get_libcpptraj_files(self):
        """Find libcpptraj library files"""
        libcpptraj = self.platform.get_library_extension()

        # Try OpenMP version first
        if self.config.openmp_flag:
            libcpptraj_files = glob(os.path.join(self.cpptraj_info.lib_dir, 'libcpptraj_omp') + '*')
            if libcpptraj_files:
                return libcpptraj_files

        # Fall back to regular version
        return glob(os.path.join(self.cpptraj_info.lib_dir, libcpptraj) + '*')

    def get_extensions(self):
        """Build and return extension modules"""
        if not self.config.compile_extensions or self.config.is_sdist:
            return []

        # Set environment variables for compilers
        setenv_cc_cxx(self.cpptraj_info.ambertools_distro, [], [])

        # Get compile arguments
        extra_compile_args, extra_link_args = self.platform.get_compile_args(self.config.debug)

        # Add platform-specific rpath
        self.platform.add_rpath_if_linux(
            self.cpptraj_info.lib_dir, extra_link_args, extra_compile_args
        )

        return get_ext_modules(
            cpptraj_info=self.cpptraj_info,
            pytraj_src=self.pytraj_src,
            compile_c_extension=self.config.compile_extensions,
            is_released=False,  # Based on original logic
            need_cython=self.need_cython,
            cpptraj_included=self.cpptraj_included,
            libcpptraj_files=self.libcpptraj_files,
            openmp_flag=self.config.openmp_flag,
            use_amberlib=self.config.use_amberlib,
            cython_directives=self.config.get_cython_directives(),
            Extension=Extension,
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args,
            define_macros=self.config.get_define_macros(),
            use_pip=self.config.use_pip,
            tarfile=self.config.is_sdist,
            use_prebuilt_cythonized_files=self.config.use_prebuilt
        )

    def get_package_data(self):
        """Get package data files"""
        return get_package_data(use_pip=self.config.use_pip)