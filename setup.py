from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources

from beastgenpy import __version__, _program

setup(name='ABSynthE',
      version=__version__,
      packages=find_packages(),
      scripts=[
            ],
      install_requires=[
        ],
      description='Agent based synthetic epidemic',
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
