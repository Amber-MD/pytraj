#!/usr/bin/env python
"""
Refactored setup.py - cleaner structure focusing on main setup logic
Complex build logic moved to separate modules for better maintainability
"""

import os
import sys
from setuptools import setup, Extension

# Import build configuration and utilities
from build_config import BuildConfig
from platform_utils import PlatformHandler
from extension_builder import ExtensionBuilder

def main():
    """Main setup entry point with cleaner flow"""

    # Parse command line arguments and flags
    config = BuildConfig()

    # Handle platform-specific setup
    platform = PlatformHandler()
    platform.configure_environment(config)

    # Build extensions if needed
    ext_builder = ExtensionBuilder(config, platform)
    ext_modules = ext_builder.get_extensions() if config.compile_extensions else []

    # Get package data
    package_data = ext_builder.get_package_data()

    setup(
        name="pytraj",
        version=config.version,
        author="Hai Nguyen",
        url="https://github.com/Amber-MD/pytraj",
        packages=config.packages,
        description="Python API for cpptraj: a data analysis package for biomolecular simulation",
        license="GPL v3",
        install_requires=['numpy'],
        classifiers=config.classifiers,
        ext_modules=ext_modules,
        package_data={'pytraj': package_data},
        cmdclass=config.cmdclass,
    )

if __name__ == "__main__":
    main()