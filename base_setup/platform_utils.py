"""
Platform-specific configuration and environment setup
"""
import os
import sys
import subprocess
import warnings

# Import from local base_setup module
from .base_setup import is_clang

class PlatformHandler:
    """Handle platform-specific build configuration"""

    def __init__(self):
        self.is_darwin = sys.platform.startswith('darwin')
        self.is_linux = sys.platform.startswith('linux')
        self.is_windows = sys.platform.startswith('win')

    def configure_environment(self, config):
        """Configure build environment based on platform"""
        if self.is_darwin:
            self._configure_macos(config)
        elif self.is_linux:
            self._configure_linux(config)
        elif self.is_windows:
            self._configure_windows(config)

        # Platform-specific OpenMP handling
        if self.is_darwin or self.is_windows:
            print(f'pytraj does not support openmp on {sys.platform} - disabled')
            config.openmp_flag = ''
            config.disable_openmp = True

    def _configure_macos(self, config):
        """macOS-specific configuration"""
        # Set default compilers
        if not os.getenv('CXX'):
            os.environ['CXX'] = 'clang++'
        if not os.getenv('CC'):
            os.environ['CC'] = 'clang'

        # Set deployment target if using clang
        if is_clang(os.getenv('CXX', 'clang++')):
            self._set_macos_deployment_target()

    def _configure_linux(self, config):
        """Linux-specific configuration"""
        # Set default compilers if not set
        if not os.getenv('CC'):
            os.environ['CC'] = 'gcc'
        if not os.getenv('CXX'):
            os.environ['CXX'] = 'g++'

    def _configure_windows(self, config):
        """Windows-specific configuration"""
        # Windows-specific compiler settings would go here
        pass

    def _set_macos_deployment_target(self):
        """Set appropriate macOS deployment target"""
        darwin_to_osx = {
            '11': '10.7', '12': '10.8', '13': '10.9', '14': '10.10',
            '15': '10.11', '16': '10.12', '17': '10.13', '18': '10.14',
            '19': '10.15', '20': '11.5', '21': '12.2', '22': '13.0', '23': '14.0'
        }

        darwin_major = os.uname()[2].split('.')[0]
        if darwin_major in darwin_to_osx:
            os.environ['MACOSX_DEPLOYMENT_TARGET'] = darwin_to_osx[darwin_major]
        else:
            warnings.warn("darwin version mapping needs update")
            try:
                result = subprocess.run(
                    ["sw_vers", "-productVersion"],
                    capture_output=True, text=True
                )
                version = ".".join(result.stdout.strip().split(".")[:2])
                os.environ["MACOSX_DEPLOYMENT_TARGET"] = version
            except Exception:
                pass

    def get_library_extension(self):
        """Get platform-appropriate library extension"""
        if self.is_darwin:
            return 'libcpptraj.dylib'
        elif self.is_windows:
            return 'libcpptraj.dll'
        else:
            return 'libcpptraj.so'

    def get_compile_args(self, debug=False):
        """Get platform-specific compile arguments"""
        if self.is_windows:
            return [], []  # No special args for Windows
        else:
            if debug:
                return ['-O0', '-ggdb'], ['-O0', '-ggdb']
            else:
                return [], []

    def add_rpath_if_linux(self, lib_dir, extra_link_args, extra_compile_args):
        """Add rpath for Linux builds"""
        if self.is_linux:
            rpath_arg = f'-Wl,-rpath={lib_dir}'
            extra_link_args.append(rpath_arg)
            extra_compile_args.append(rpath_arg)
            print(f'set rpath to {lib_dir}')

    def get_libraries(self):
        """Get platform-appropriate library names"""
        if self.is_windows:
            return ['libcpptraj']
        else:
            return ['cpptraj']

    def get_additional_include_dirs(self):
        """Get additional include directories for platform"""
        if self.is_windows:
            # Add unistd.h for Windows
            return [os.path.join(os.path.dirname(__file__), '..', 'scripts', '3rd_party')]
        return []