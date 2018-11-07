"""
List of test cases for NetworkEX

Kiri Choi 
"""

import os

dir_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'testcases')
testfiles = [f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))]

content = []
for i in testfiles:
    sbmlstr = open(os.path.join(dir_path, i), 'r')
    content.append(sbmlstr.read())
    sbmlstr.close()

list_models = [s for s in testfiles if "test_list" in s]
_LIST = []
for i in range(len(list_models)):
    sbmlstr = open(os.path.join(dir_path, list_models[i]), 'r')
    _LIST.append(sbmlstr.read())
    sbmlstr.close()

class testmodels():

    BIOMD0000000010 = content[testfiles.index('test_BIOMD0000000010.xml')]
    BIBI = content[testfiles.index('test_bibi.xml')]
    BRANCHED = content[testfiles.index('test_branched.xml')]
    REVERSIBLE = content[testfiles.index('test_reversible.xml')]
    SYMPY = content[testfiles.index('test_sympy.xml')]
    FEEDBACK = content[testfiles.index('test_feedback.xml')]
    ASPARTATE = content[testfiles.index('test_AspartateMetabolism.xml')]
    LIST = _LIST
    
        
        