# -*- coding: utf-8 -*-
# bases/algorithm/algorithmbase_test.py

from __future__ import absolute_import # PEP328
from builtins import zip
from nose.tools import raises, assert_raises
from numpy import array as np_array
from numpy import uint32, float64, dtype
from bases.algorithm.algorithmbase import (AlgorithmNameError,
                                           AlgorithmParameterError)
from bases.algorithm import (AlgorithmBase, Parameter)

TestPar = Parameter("testPar", 5.0, valueRange = (4.0, 13.0))

@raises(AlgorithmNameError)
def testName():
    """name is mandatory"""
    class DummyAlgo(AlgorithmBase):
        pass
    t = DummyAlgo()

def testParam():
    """Algorithm without parameters allowed"""
    class DummyAlgo(AlgorithmBase):
        pass
    atype = DummyAlgo.factory("testalgo")
    t = atype()

def testParam1():
    class DummyAlgo(AlgorithmBase):
        pass
    atype = DummyAlgo.factory("testalgo", TestPar)
    t = atype()

@raises(AlgorithmParameterError)
def testParam4():
    """Provided parameters have to be a subclass of ParameterBase"""
    class DummyAlgo(AlgorithmBase):
        pass
    atype = DummyAlgo.factory("testalgo", "dummy")
    t = atype()

def testTypeVsInstance():
    class DummyAlgo(AlgorithmBase):
        pass
    atype = DummyAlgo.factory("testalgo", TestPar)
    ainst = atype()
    # same name
    assert atype.name() == ainst.name()
    # same (default) value
    assert atype.testPar.value() == ainst.testPar.value()
    # class vs instance
    assert type(atype.testPar) is type
    assert type(ainst.testPar) is TestPar
    assert id(atype.testPar) != id(ainst.testPar)
    # type of instance should remain in class too
    assert id(atype.testPar) == id(type(ainst).testPar)

class DummyAlgo(AlgorithmBase):
    shortName = "Dummy"
    parameters = (TestPar, Parameter(name = "test", value = 3.4,
                                     valueRange = (1,5)))

    def __init__(self):
        super(DummyAlgo, self).__init__()
        self.test.setOnValueUpdate(self.dummy)

    def dummy(self, v):
        pass

def testSerialize():
    import pickle
    DummyAlgo.factory()
    da = DummyAlgo()
    data = pickle.dumps(da)
    db = pickle.loads(data)
    assert da == db
    assert id(da) != id(db)
    # make sure the parameter instances are the same
    # in different places of the class
    for p0, p1 in zip(db.params(), (db.testPar, db.test)):
        assert id(p0) == id(p1), \
                "Ids do not match for '{}'!".format(p0.name())

# vim: set ts=4 sts=4 sw=4 tw=0:
