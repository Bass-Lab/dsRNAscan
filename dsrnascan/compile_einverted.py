#!/usr/bin/env python3
"""
Compile einverted with G-U wobble patch if needed.
This can be called directly or during installation.
"""

import os
import sys
import subprocess
import shutil
import platform

def get_platform_binary_name():
    """Get the expected binary name for the current platform"""
    system = platform.system().lower()
    machine = platform.machine().lower()
    
    if system == 'darwin':
        if 'arm' in machine or 'aarch64' in machine:
            return 'einverted_darwin_arm64'
        else:
            return 'einverted_darwin_x86_64'
    elif system == 'linux':
        if 'arm' in machine or 'aarch64' in machine:
            return 'einverted_linux_aarch64'
        else:
            return 'einverted_linux_x86_64'
    elif system == 'windows':
        return 'einverted_windows_x86_64.exe'
    else:
        return 'einverted'

def find_tools_dir():
    """Find the tools directory"""
    # Try different possible locations
    possible_paths = [
        os.path.join(os.path.dirname(__file__), 'tools'),
        os.path.join(os.path.dirname(os.path.dirname(__file__)), 'dsrnascan', 'tools'),
        os.path.join(os.getcwd(), 'dsrnascan', 'tools'),
    ]
    
    for path in possible_paths:
        if os.path.exists(path):
            return path
    
    # Create it if it doesn't exist
    default_path = os.path.join(os.path.dirname(__file__), 'tools')
    os.makedirs(default_path, exist_ok=True)
    return default_path

def compile_einverted():
    """Compile einverted with G-U patch"""
    tools_dir = find_tools_dir()
    target_binary = os.path.join(tools_dir, 'einverted')
    
    # Check if already exists
    if os.path.exists(target_binary):
        print(f"einverted already exists at {target_binary}")
        return True
    
    # Check for platform-specific binary
    platform_binary_name = get_platform_binary_name()
    platform_binary = os.path.join(tools_dir, platform_binary_name)
    
    if os.path.exists(platform_binary):
        print(f"Using pre-compiled binary: {platform_binary_name}")
        shutil.copy2(platform_binary, target_binary)
        os.chmod(target_binary, 0o755)
        return True
    
    # Try to compile
    print(f"No pre-compiled binary found for {platform.system()} {platform.machine()}")
    print("Attempting to compile einverted with G-U wobble patch...")
    
    # Find compilation files
    package_root = os.path.dirname(os.path.dirname(__file__))
    compile_script = os.path.join(package_root, 'compile_patched_einverted.sh')
    patch_file = os.path.join(package_root, 'einverted.patch')
    
    if not os.path.exists(compile_script) or not os.path.exists(patch_file):
        print(f"ERROR: Cannot find compilation files")
        print(f"  Script: {compile_script} (exists: {os.path.exists(compile_script)})")
        print(f"  Patch: {patch_file} (exists: {os.path.exists(patch_file)})")
        print("\neinverted binary is required for dsRNAscan to work properly.")
        print("Please compile it manually or install on a supported platform.")
        return False
    
    try:
        # Make script executable
        os.chmod(compile_script, 0o755)
        
        # Run compilation
        result = subprocess.run(
            ['bash', compile_script],
            cwd=package_root,
            capture_output=True,
            text=True
        )
        
        if result.returncode == 0:
            print("Successfully compiled einverted with G-U wobble patch")
            return True
        else:
            print(f"Compilation failed: {result.stderr}")
            return False
            
    except Exception as e:
        print(f"Error during compilation: {e}")
        return False

if __name__ == "__main__":
    success = compile_einverted()
    sys.exit(0 if success else 1)