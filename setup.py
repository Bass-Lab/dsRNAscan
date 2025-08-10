import os
import sys
import subprocess
import platform
from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools.command.develop import develop

# Version information
__version__ = '0.3.0'

class CustomInstallCommand(install):
    """Custom installation to handle einverted binary"""
    def run(self):
        self.setup_einverted()
        install.run(self)
        
    def setup_einverted(self):
        """Set up platform-specific einverted binary"""
        tools_dir = os.path.join(os.path.dirname(__file__), 'tools')
        
        # Detect platform
        system = platform.system().lower()
        machine = platform.machine().lower()
        
        # Map platform to binary name
        if system == 'darwin':
            if 'arm' in machine or 'aarch64' in machine:
                binary_name = 'einverted_darwin_arm64'
            else:
                binary_name = 'einverted_darwin_x86_64'
        elif system == 'linux':
            if 'arm' in machine or 'aarch64' in machine:
                binary_name = 'einverted_linux_aarch64'
            else:
                binary_name = 'einverted_linux_x86_64'
        elif system == 'windows':
            binary_name = 'einverted_windows_x86_64.exe'
        else:
            binary_name = None
            
        # Check for platform-specific binary
        if binary_name:
            platform_binary = os.path.join(tools_dir, binary_name)
            target_binary = os.path.join(tools_dir, 'einverted')
            
            if os.path.exists(platform_binary):
                # Copy platform-specific binary to einverted
                import shutil
                shutil.copy2(platform_binary, target_binary)
                os.chmod(target_binary, 0o755)
                print(f"Using precompiled binary: {binary_name}")
                return
                
        # Fallback: check if generic einverted exists
        einverted_path = os.path.join(tools_dir, 'einverted')
        if os.path.exists(einverted_path):
            print(f"Using existing einverted binary")
            os.chmod(einverted_path, 0o755)
        else:
            print("Warning: No suitable einverted binary found for your platform")
            print(f"Platform: {system} {machine}")
            print("Please compile einverted from source or use conda installation")

class CustomDevelopCommand(develop):
    """Custom develop command to handle einverted binary"""
    def run(self):
        # Use the same setup method as install
        installer = CustomInstallCommand(self.distribution)
        installer.setup_einverted()
        develop.run(self)

# Read long description from README
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='dsrnascan',
    version=__version__,
    author='Bass Lab',
    author_email='',
    description='A tool for genome-wide prediction of double-stranded RNA structures',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/Bass-Lab/dsRNAscan',
    project_urls={
        "Bug Tracker": "https://github.com/Bass-Lab/dsRNAscan/issues",
        "Documentation": "https://github.com/Bass-Lab/dsRNAscan/blob/main/README.md",
        "Source Code": "https://github.com/Bass-Lab/dsRNAscan",
    },
    py_modules=['dsRNAscan', 'dsRNAscan-mpi'],
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
    ],
    python_requires='>=3.7',
    install_requires=[
        'biopython>=1.78',
        'numpy>=1.19',
        'pandas>=1.1',
        'ViennaRNA>=2.4',
    ],
    extras_require={
        'mpi': ['mpi4py>=3.0', 'parasail>=1.2'],
        'dev': [
            'pytest>=6.0',
            'pytest-cov>=2.0',
            'black>=22.0',
            'flake8>=4.0',
            'mypy>=0.900',
        ],
    },
    entry_points={
        'console_scripts': [
            'dsrnascan=dsRNAscan:main',
            # 'dsrnascan-mpi=dsRNAscan_mpi:main',  # Uncomment when module is renamed
        ],
    },
    include_package_data=True,
    package_data={
        '': ['tools/einverted', 'einverted.patch', 'einverted.c'],
    },
    cmdclass={
        'install': CustomInstallCommand,
        'develop': CustomDevelopCommand,
    },
    zip_safe=False,
)