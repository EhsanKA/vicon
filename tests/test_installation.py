#!/usr/bin/env python3

import sys
import os
import importlib
import subprocess
from pathlib import Path

def check_import(module_name):
    """Check if a module can be imported."""
    try:
        importlib.import_module(module_name)
        print(f"✓ Successfully imported {module_name}")
        return True
    except ImportError as e:
        print(f"✗ Failed to import {module_name}: {str(e)}")
        return False

def check_command(command):
    """Check if a command is available and executable."""
    try:
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"✓ Command '{command[0]}' is available")
            return True
        else:
            print(f"✗ Command '{command[0]}' failed: {result.stderr}")
            return False
    except FileNotFoundError:
        print(f"✗ Command '{command[0]}' not found")
        return False

def main():
    print("Testing VICON installation...\n")
    
    # Check Python version
    print(f"Python version: {sys.version}")
    
    # Check VICON modules
    modules = [
        "vicon",
        "vicon.dereplication",
        "vicon.alignment",
        "vicon.processing",
        "vicon.visualization",
        "vicon.io",
        "vicon.utils"
    ]
    
    print("\nChecking VICON modules:")
    module_results = [check_import(module) for module in modules]
    
    # Check required commands
    print("\nChecking required commands:")
    commands = [
        ["vicon-run", "--help"],
        ["viralmsa", "--help"],
        ["vsearch", "--version"]
    ]
    command_results = [check_command(cmd) for cmd in commands]
    
    # Summary
    print("\nInstallation Test Summary:")
    print(f"VICON modules: {sum(module_results)}/{len(module_results)} passed")
    print(f"Commands: {sum(command_results)}/{len(command_results)} passed")
    
    if all(module_results) and all(command_results):
        print("\n✓ VICON installation appears to be working correctly!")
    else:
        print("\n✗ Some tests failed. Please check the output above for details.")

if __name__ == "__main__":
    main() 