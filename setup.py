import os
import sys
import subprocess
import platform
from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools.command.build_py import build_py

# Version information
__version__ = '0.3.0'

class CustomBuildPy(build_py):
    """Custom build command to compile einverted during build phase"""
    def run(self):
        # Compile einverted before building the package
        self.setup_einverted()
        build_py.run(self)
    
    def setup_einverted(self):
        """Set up platform-specific einverted binary with G-U wobble patch"""
        import shutil
        import tempfile
        import urllib.request
        
        tools_dir = os.path.join(os.path.dirname(__file__), 'dsrnascan', 'tools')
        os.makedirs(tools_dir, exist_ok=True)
        
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
            binary_name = 'einverted'
            
        # First check if we have a pre-compiled patched binary
        platform_binary = os.path.join(tools_dir, binary_name)
        target_binary = os.path.join(tools_dir, 'einverted')
        
        if os.path.exists(platform_binary):
            # Check if it's actually a binary (not a script)
            import stat
            with open(platform_binary, 'rb') as f:
                header = f.read(4)
            # Check for binary headers (ELF, Mach-O, PE)
            is_binary = header in [
                b'\x7fELF',  # Linux ELF
                b'\xcf\xfa\xed\xfe',  # macOS Mach-O 64-bit
                b'\xce\xfa\xed\xfe',  # macOS Mach-O 32-bit
                b'\xca\xfe\xba\xbe',  # macOS Universal binary
                b'MZ\x90\x00',  # Windows PE
            ]
            
            if is_binary:
                shutil.copy2(platform_binary, target_binary)
                os.chmod(target_binary, 0o755)
                print(f"Using precompiled patched binary: {binary_name}")
                return
            else:
                print(f"Found {binary_name} but it's not a valid binary, will compile from source")
                
        # First try minimal placeholder compilation
        print("No pre-compiled binary found. Setting up einverted...")
        
        original_dir = os.getcwd()
        
        # Try to compile a minimal placeholder first
        try:
            minimal_c = os.path.join(os.path.dirname(__file__), 'compile_minimal_einverted.c')
            if os.path.exists(minimal_c):
                print("Compiling minimal einverted placeholder...")
                subprocess.run(['gcc', '-o', target_binary, minimal_c], 
                             check=True, stderr=subprocess.DEVNULL)
                os.chmod(target_binary, 0o755)
                print("✓ Installed einverted placeholder. For full functionality, install EMBOSS.")
                return
        except (subprocess.CalledProcessError, FileNotFoundError):
            pass  # gcc not available, try other methods
        
        # Try to compile from EMBOSS source with patch
        compile_script = os.path.join(os.path.dirname(__file__), 'compile_patched_einverted.sh')
        patch_file = os.path.join(os.path.dirname(__file__), 'einverted.patch')
        
        # Try to compile if we have the compilation script or user explicitly requests it
        # Always try to compile during pip install to ensure G-U patch is applied
        if os.path.exists(compile_script) and os.path.exists(patch_file):
            print("Attempting to compile einverted with G-U wobble patch...")
            try:
                # Make script executable
                os.chmod(compile_script, 0o755)
                
                # Run the compilation script
                result = subprocess.run(
                    ['bash', compile_script],
                    cwd=os.path.dirname(__file__),
                    capture_output=True,
                    text=True
                )
                
                if result.returncode == 0:
                    print("✓ Successfully compiled einverted with G-U wobble patch")
                    
                    # Check if the binary was created
                    if os.path.exists(target_binary):
                        os.chmod(target_binary, 0o755)
                        return
                    else:
                        print("WARNING: Compilation succeeded but binary not found")
                else:
                    print(f"WARNING: Compilation failed with exit code {result.returncode}")
                    if result.stderr:
                        print(f"Error output: {result.stderr}")
                        
            except Exception as e:
                print(f"WARNING: Could not compile einverted with patch: {e}")
                print("This may be due to missing dependencies (gcc, make, etc.)")
                print("You can try installing EMBOSS separately or compiling manually")
            
        # Last resort: Create a placeholder binary if needed
        if not os.path.exists(target_binary):
            print("WARNING: Creating placeholder einverted binary")
            print("Note: einverted will need to be installed separately for full functionality")
            # Create a simple script that explains the issue
            with open(target_binary, 'w') as f:
                f.write("#!/bin/sh\n")
                f.write('echo "Error: einverted binary not properly installed."\n')
                f.write('echo "Please install EMBOSS or use conda: conda install -c bioconda emboss"\n')
                f.write('exit 1\n')
            os.chmod(target_binary, 0o755)

class CustomInstallCommand(install):
    """Custom installation to use CustomBuildPy"""
    def run(self):
        install.run(self)

class CustomDevelopCommand(develop):
    """Custom develop command to handle einverted binary"""
    def run(self):
        # Use the build command to compile einverted
        self.run_command('build_py')
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
    packages=['dsrnascan'],
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
            'dsrnascan=dsrnascan:main',
        ],
    },
    include_package_data=True,
    package_data={
        'dsrnascan': ['tools/*'],
        '': ['einverted.patch', 'compile_patched_einverted.sh', 'compile_minimal_einverted.c'],
    },
    cmdclass={
        'build_py': CustomBuildPy,
        'install': CustomInstallCommand,
        'develop': CustomDevelopCommand,
    },
    zip_safe=False,
)