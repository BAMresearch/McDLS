# -*- coding: utf-8 -*-
# algorithmbase_test.py

from algorithmbase import AlgorithmBase, AlgorithmNameError, AlgorithmParameterError
from parameter import ParameterNumerical, ParameterError
from nose.tools import raises, assert_raises
from numpy import array as np_array
from numpy import uint32, float64, dtype

class TestParam(ParameterNumerical):
    name = "testPar"
    defaultValue = 5.0
    valueRange = (4.0, 13.0)

def testParam1():
    class TestAlgo(AlgorithmBase):
        shortName = "testalgo"
    t = TestAlgo()

def testParam2():
    class TestAlgo(AlgorithmBase):
        shortName = "testalgo"
        parameters = TestParam
    t = TestAlgo()

@raises(AlgorithmNameError)
def testName():
    class TestAlgo(AlgorithmBase):
        pass
    t = TestAlgo()

@raises(AlgorithmParameterError)
def testParam3():
    AlgorithmBase.parameters = TestParam
    class TestAlgo(AlgorithmBase):
        shortName = "testalgo"
    t = TestAlgo()

@raises(AlgorithmParameterError)
def testParam4():
    AlgorithmBase.parameters = None
    class TestAlgo(AlgorithmBase):
        shortName = "testalgo"
        parameters = "test"
    t = TestAlgo()

@raises(AlgorithmParameterError)
def testParam5():
    AlgorithmBase.parameters = None
    class TestAlgo(AlgorithmBase):
        shortName = "testalgo"
        parameters = TestParam
        testPar = None
    t = TestAlgo()

# vim: set ts=4 sts=4 sw=4 tw=0:
