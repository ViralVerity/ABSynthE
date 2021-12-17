from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources

from absynthe import __version__, _program

setup(name='absynthe',
      version=__version__,
      packages=find_packages(),
      description='Agent based synthetic epidemic',
      install_requires=[
            "numpy>=1.19.4",
            "scipy>=1.4.1"
            
      ],
      url='https://github.com/ViralVerity/ABSynthE',
      author='Verity Hill',
      author_email='verity.hill@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = absynthe.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
