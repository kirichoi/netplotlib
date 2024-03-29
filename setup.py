# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 21:48:16 2018

@author: kirichoi
"""

from setuptools import setup
import os

with open(os.path.join(os.path.dirname(__file__), 'netplotlib/VERSION.txt'), 'r') as f:
    _version = f.read().rstrip()

setup(name='netplotlib',
      version=_version,
      author='Kiri Choi',
      description='netplotlib: A simple package for visualizing reaction network models',
      long_description='netplotlib is a simple extension to NetworkX and matplotlib to draw reaction network diagrams from SBML or Antimony strings with ease. netplotlib also supports drawing weighted reaction network diagrams from an ensemble of models based on reaction frequency or user-supplied weights.',
      packages=[
          'netplotlib',
		  'netplotlib.testcases',
      ],
      package_data={
          'netplotlib.testcases': ['*.xml'],
		  'netplotlib': ['*.txt'],
      },
      install_requires=[
        'tellurium>=2.1.0',
        'networkx>=2.1',
        'numpy>=1.15.0',
        'scipy>=0.19.0',
        'matplotlib>=2.0.2',
        'libroadrunner>=1.4.24',
        'antimony>=2.9.4',
        'sympy>=1.1.1',
        'python-libsbml>=5.15.0',
        'opencv-python>=4.5.1.48',
        ],
      classifiers=[
        'License :: OSI Approved :: MIT License',
        'Development Status :: 5 - Production/Stable',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'
        ],
      )
