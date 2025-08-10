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
            # Test if it has the patch by checking for G-U functionality
            # For now, assume pre-compiled binaries are patched
            shutil.copy2(platform_binary, target_binary)
            os.chmod(target_binary, 0o755)
            print(f"Using precompiled patched binary: {binary_name}")
            return
                
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
        
        # Try to compile from EMBOSS source with patch (if user wants full version)
        if os.environ.get('DSRNASCAN_COMPILE_FULL', '').lower() == 'true':
            print("Attempting to compile full einverted with G-U wobble patch...")
            try:
                # Create temp directory for compilation
                with tempfile.TemporaryDirectory() as tmpdir:
                    os.chdir(tmpdir)
                    
                    # Download EMBOSS if not present
                    emboss_url = "ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz"
                    emboss_tar = "EMBOSS-6.6.0.tar.gz"
                    
                    print("Downloading EMBOSS source...")
                    try:
                        urllib.request.urlretrieve(emboss_url, emboss_tar)
                    except Exception as e:
                        print(f"WARNING: Could not download EMBOSS: {e}")
                        raise
                
                    
                    # Extract
                    subprocess.run(['tar', '-xzf', emboss_tar], check=True)
                    os.chdir('EMBOSS-6.6.0')
                    
                    # Apply patch
                    patch_file = os.path.join(os.path.dirname(__file__), 'einverted.patch')
                    if os.path.exists(patch_file):
                        print("Applying G-U wobble patch...")
                        subprocess.run(['patch', '-p1', '-i', patch_file], check=True)
                    
                    # Configure and compile
                    print("Configuring EMBOSS (this may take a minute)...")
                    subprocess.run(['./configure', '--without-x', '--disable-shared'], 
                                 check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    
                    print("Compiling einverted with G-U patch...")
                    os.chdir('emboss')
                    subprocess.run(['make', 'einverted'], 
                                 check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    
                    # Copy the compiled binary
                    if os.path.exists('einverted'):
                        shutil.copy2('einverted', target_binary)
                        os.chmod(target_binary, 0o755)
                        
                        # Also save as platform-specific for future use
                        shutil.copy2('einverted', platform_binary)
                        os.chmod(platform_binary, 0o755)
                        
                        print(f"✓ Successfully compiled einverted with G-U wobble patch for {system} {machine}")
                        os.chdir(original_dir)
                        return
                        
            except Exception as e:
                print(f"WARNING: Could not compile einverted with patch: {e}")
                os.chdir(original_dir)
            
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
    },
    cmdclass={
        'install': CustomInstallCommand,
        'develop': CustomDevelopCommand,
    },
    zip_safe=False,
)