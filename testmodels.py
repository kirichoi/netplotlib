"""
List of test models for NetworkEX

Kiri Choi 
"""

import os

dir_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'testcases')
testfiles = [f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))]


content = []
for i in testfiles:
    sbmlstr = open(os.path.join(dir_path, i), 'r')
    content.append(sbmlstr.read())

class testmodels:
    
    BIOMD0000000010 = content[0]
    BIBI = content[1]
    BRANCHED = content[2]
    LIST = content[3:6]
    REVERSIBLE = content[6]
    