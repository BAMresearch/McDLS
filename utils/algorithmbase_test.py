# -*- coding: utf-8 -*-
# algorithmbase_test.py

from algorithmbase import AlgorithmBase
from parameter import ParameterNumerical, ParameterError
from nose.tools import raises, assert_raises
from numpy import array as np_array
from numpy import uint32, float64, dtype

class TestAlgo(AlgorithmBase):
    shortName = "filtertest"
    parameters = [ParameterNumerical("testPar", 5.0, [4.0, 13.0]),]
    def calc(self, data):
        pass

def testCopy():
    tf = TestAlgo()
    cf = tf.copy()
    assert id(cf) != id(tf) # different instances
    assert cf == tf
    cf.testPar = 6.0 # change arbitrarily
    assert cf != tf
    tf.testPar = cf.testPar
    cf.testPar = 5.0 # reset to default value
    assert cf != tf

# vim: set ts=4 sts=4 sw=4 tw=0:
