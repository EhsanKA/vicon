import os
import platform
import shutil
from setuptools import setup, find_packages
from setuptools.command.install import install

class PostInstallCommand(install):
    """Custom post-installation script to download vsearch based on architecture."""
    def run(self):
        install.run(self)
        arch = platform.machine()
        if arch == "x86_64":
            url = "https://github.com/torognes/vsearch/releases/download/v2.30.0/vsearch-2.30.0-linux-x86_64"
        else:
            raise RuntimeError(f"Unsupported architecture: {arch}")
        
        bin_dir = os.path.join(self.install_scripts, "..", "bin")
        vsearch_dest = os.path.join(bin_dir, "vsearch")
        os.system(f"wget {url} -O {vsearch_dest}")
        os.chmod(vsearch_dest, 0o755)

setup(
    name="vicon",
    version="1.0.3",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy>=2.2.6",
        "pandas>=2.2.3",
        "matplotlib>=3.10.3",
        "pyyaml>=6.0.2",
        "biopython>=1.85",
        "setuptools",
        "pip",
    ],
    entry_points={
        "console_scripts": [
            "viralmsa=vicon.scripts.ViralMSA:main",
        ],
    },
    cmdclass={
        "install": PostInstallCommand,
    },
)