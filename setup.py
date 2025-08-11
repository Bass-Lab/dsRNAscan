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
                
        # No pre-compiled binary found
        print(f"WARNING: No pre-compiled einverted binary found for {system} {machine}")
        print("The package will install, but einverted compilation may be needed.")
        
        # Try to compile from EMBOSS source with patch
        # Use absolute path resolution that works in sdist builds
        setup_dir = os.path.dirname(os.path.abspath(__file__))
        
        # Try multiple possible locations for the files
        possible_dirs = [
            setup_dir,
            os.getcwd(),
            os.path.join(os.getcwd(), 'dsrnascan'),
            os.path.dirname(setup_dir),
        ]
        
        compile_script = None
        patch_file = None
        
        for dir_path in possible_dirs:
            test_script = os.path.join(dir_path, 'compile_patched_einverted.sh')
            test_patch = os.path.join(dir_path, 'einverted.patch')
            if os.path.exists(test_script) and os.path.exists(test_patch):
                compile_script = test_script
                patch_file = test_patch
                setup_dir = dir_path
                break
        
        if not compile_script:
            # Use default paths for error reporting
            compile_script = os.path.join(setup_dir, 'compile_patched_einverted.sh')
            patch_file = os.path.join(setup_dir, 'einverted.patch')
        
        # Try to compile if we're in a development environment (not a pip install)
        if os.path.exists(compile_script) and os.path.exists(patch_file) and os.environ.get('DSRNASCAN_COMPILE', 'false').lower() == 'true':
            print(f"Found compilation files:")
            print(f"  Script: {compile_script}")
            print(f"  Patch: {patch_file}")
            print("Compiling einverted with G-U wobble patch...")
            try:
                # Make script executable
                os.chmod(compile_script, 0o755)
                
                # Run the compilation script
                # Set environment to help the script find files
                env = os.environ.copy()
                env['SETUP_DIR'] = setup_dir
                env['TARGET_DIR'] = tools_dir
                
                result = subprocess.run(
                    ['bash', compile_script],
                    cwd=setup_dir,
                    capture_output=True,
                    text=True,
                    env=env
                )
                
                if result.returncode == 0:
                    print("✓ Successfully compiled einverted with G-U wobble patch")
                    
                    # Check if the binary was created
                    if os.path.exists(target_binary):
                        os.chmod(target_binary, 0o755)
                        return
                    else:
                        raise RuntimeError("Compilation succeeded but binary not found at " + target_binary)
                else:
                    error_msg = f"Compilation failed with exit code {result.returncode}"
                    if result.stderr:
                        error_msg += f"\nError output: {result.stderr}"
                    if result.stdout:
                        error_msg += f"\nOutput: {result.stdout}"
                    raise RuntimeError(error_msg)
                        
            except Exception as e:
                print(f"ERROR: Could not compile einverted with G-U patch: {e}")
                print("\n" + "="*60)
                print("IMPORTANT: dsRNAscan requires einverted with G-U wobble patch")
                print("Please ensure you have the following installed:")
                print("  - gcc/clang compiler")
                print("  - make")
                print("  - patch")
                print("  - wget or curl")
                print("\nOr manually compile by running:")
                print("  ./compile_patched_einverted.sh")
                print("="*60 + "\n")
                raise RuntimeError(f"Failed to compile einverted with G-U patch: {e}")
        else:
            # Don't fail the build, just warn
            print(f"INFO: einverted will need to be compiled after installation")
            print(f"  Run: python -m dsrnascan.compile_einverted")
            print(f"  Or set DSRNASCAN_COMPILE=true during installation to compile automatically")
            
            # Create a placeholder script if no binary exists
            if not os.path.exists(target_binary):
                with open(target_binary, 'w') as f:
                    f.write('#!/bin/sh\n')
                    f.write('echo "ERROR: einverted binary not compiled yet."\n')
                    f.write('echo "Please run: python -m dsrnascan.compile_einverted"\n')
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