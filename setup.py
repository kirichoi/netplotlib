# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 21:48:16 2018

@author: kirichoi
"""

from setuptools import setup
import os

try:
    with open(os.path.join(os.path.dirname(__file__), 'VERSION.txt'), 'r') as f:
        version = f.read().rstrip()
except IOError:
    with open(os.path.join(os.path.dirname(__file__), 'networkex/VERSION.txt'), 'r') as f:
        version = f.read().rstrip()

setup(name='networkex',
      version=version,
      author='Kiri Choi',
      description='NetworkEx: A simple extension for NetworkX for Antimony models',
      packages=[
          'networkex'
      ],
      install_requires=[
        'tellurium>=2.1.0',
        'networkx>=2.1',
        'numpy>=1.15.0',
        'scipy>=0.19.0',
        'matplotlib>=2.0.2',
        'libroadrunner>=1.4.24',
        'antimony>=2.9.4',
        'sympy>=1.1.1'
        ]
      )
