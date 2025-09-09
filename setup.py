import os
import sys
import subprocess
import platform
from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools.command.build_py import build_py

# Version information
__version__ = '0.4.1'

def compile_einverted_now():
    """Compile einverted immediately when setup.py runs"""
    print("\n" + "="*60)
    print("dsRNAscan setup: Checking for einverted...")
    print("="*60)
    
    # Determine paths
    setup_dir = os.path.dirname(os.path.abspath(__file__))
    tools_dir = os.path.join(setup_dir, 'dsrnascan', 'tools')
    os.makedirs(tools_dir, exist_ok=True)
    
    target_binary = os.path.join(tools_dir, 'einverted')
    
    # Check if already compiled
    if os.path.exists(target_binary):
        with open(target_binary, 'rb') as f:
            header = f.read(4)
        if header in [b'\x7fELF', b'\xcf\xfa', b'\xce\xfa', b'\xca\xfe']:
            size = os.path.getsize(target_binary)
            print(f"✓ Found existing einverted binary ({size} bytes)")
            return
    
    # Need to compile
    compile_script = os.path.join(setup_dir, 'compile_patched_einverted.sh')
    patch_file = os.path.join(setup_dir, 'einverted.patch')
    
    if not os.path.exists(compile_script) or not os.path.exists(patch_file):
        print("WARNING: Cannot find compilation files")
        print("  einverted will not be available")
        return
    
    print("Compiling einverted with G-U wobble patch...")
    print("This may take a few minutes...")
    
    try:
        # Make script executable
        os.chmod(compile_script, 0o755)
        
        # Run compilation
        env = os.environ.copy()
        env['TARGET_DIR'] = tools_dir
        
        result = subprocess.run(
            ['bash', compile_script],
            cwd=setup_dir,
            capture_output=True,
            text=True,
            env=env,
            timeout=300
        )
        
        if result.returncode == 0 and os.path.exists(target_binary):
            size = os.path.getsize(target_binary)
            print(f"✓ Successfully compiled einverted ({size} bytes)")
        else:
            print("WARNING: Compilation failed")
            print("  Install will continue but einverted won't work")
            if result.stderr:
                print(f"  Error: {result.stderr[:200]}")
                
    except Exception as e:
        print(f"WARNING: Could not compile einverted: {e}")
        print("  Install will continue but einverted won't work")

# Compile immediately when setup.py is run for install/build
if any(cmd in sys.argv for cmd in ['install', 'build', 'bdist_wheel', 'develop']):
    compile_einverted_now()

# Read long description
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
    package_data={
        'dsrnascan': ['tools/*'],
    },
    include_package_data=True,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
    ],
    python_requires='>=3.8',
    install_requires=[
        'biopython>=1.78',
        'numpy>=1.19',
        'pandas>=1.1',
        'ViennaRNA>=2.4',
        'psutil>=5.8',
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
            'dsrna-browse=dsrnascan.browse_results:main',
        ],
    },
    zip_safe=False,
)