# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 17:00:10 2018

@author: ckiri
"""
import os

try:
	with open(os.path.join(os.path.dirname(__file__), '..', 'VERSION.txt'), 'r') as f:
		version = f.read().rstrip()
except:
	with open(os.path.join(os.path.dirname(__file__), 'VERSION.txt'), 'r') as f:
		version = f.read().rstrip()

__version__ = version
