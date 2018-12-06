import os

dir_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'testcases')
testfiles = [f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))]

content = []
for i in testfiles:
    sbmlstr = open(os.path.join(dir_path, i), 'r')
    content.append(sbmlstr.read())
    sbmlstr.close()

list_models = [s for s in testfiles if "test_list1" in s]
_LIST1 = []
for i in range(len(list_models)):
    sbmlstr = open(os.path.join(dir_path, list_models[i]), 'r')
    _LIST1.append(sbmlstr.read())
    sbmlstr.close()

class testmodels():

    BIOMD0000000010 = content[testfiles.index('test_BIOMD0000000010.xml')]
    BIBI = content[testfiles.index('test_bibi.xml')]
    UNIBI = content[testfiles.index('test_unibi.xml')]
    BIUNI = content[testfiles.index('test_biuni.xml')]
    BRANCHED = content[testfiles.index('test_branched.xml')]
    REVERSIBLE = content[testfiles.index('test_reversible.xml')]
    SYMPY = content[testfiles.index('test_sympy.xml')]
    FEEDBACK = content[testfiles.index('test_feedback.xml')]
    ASPARTATE = content[testfiles.index('test_AspartateMetabolism.xml')]
    UNDEF = content[testfiles.index('test_undefBoundary.xml')]
    STOCH1 = content[testfiles.index('test_stoch1.xml')]
    STOCH2 = content[testfiles.index('test_stoch2.xml')]
    STOCH3 = content[testfiles.index('test_stoch3.xml')]
    LIST1 = _LIST1
    
        
        