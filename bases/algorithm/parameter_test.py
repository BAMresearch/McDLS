# -*- coding: utf-8 -*-
# bases/algorithm/parameter_test.py

from __future__ import absolute_import # PEP328
from bases.algorithm.parameter import (
        ParameterBase, ParameterNumerical, ParameterFloat, ParameterLog,
        factory,
        ParameterNameError, DefaultValueError, ValueRangeError,
        SuffixError, SteppingError, DecimalsError, DisplayValuesError)
from bases.algorithm import Parameter, NumberGenerator, RandomUniform
from nose.tools import raises

class Dummy(object):
    def dummyFunc(value):
        pass

def testParameterName():
    @raises(ParameterNameError)
    def testName(newName):
        p = Parameter(name, 0)
    for name in (None, "", 1.3, 0):
        yield testName, name

@raises(StandardError)
def testParameterDefaultValue1():
    p = Parameter("testpar")

@raises(DefaultValueError)
def testParameterDefaultValue2():
    p = Parameter("testpar", None)

def testParameterNumerical():
    ptype = Parameter("testpar", 3, valueRange = (1, 5),
                         suffix = "mm", stepping = 1)
    p = ptype()
    assert p.value() == 3
    assert p.valueRange() == (1, 5)
    assert p.suffix() == "mm"
    assert p.stepping() == 1
    assert p.displayValues() is None
    p.setValue(4)
    assert p.value() == 4
    assert ptype.value() == 3
    p.setValueRange((2,3))
    assert ptype.valueRange() == (1, 5)
    assert p.valueRange() == (2, 3)

def testParameterNumericalValueRange():
    @raises(ValueRangeError)
    def testValueRange(newRange):
        p = Parameter("testpar", 1, valueRange = newRange)
    for valueRange in (None, (None, 1), (1, None), (None, None), (1, 2, 3),
                       "", ("", ), (1, ""), ("", 1), ("", 1.0), (1.0, "")):
        yield testValueRange, valueRange

def testParameterNumericalSuffix():
    @raises(SuffixError)
    def testSuffix(newSuffix):
        p = Parameter("testpar", 1, valueRange = (1, 5),
                                    suffix = newSuffix)
    for suffix in ("", 1, 1.0):
        yield testSuffix, suffix

def testParameterNumericalStepping():
    @raises(SteppingError)
    def testStepping(stepping):
        p = Parameter("testpar", 1, valueRange = (1, 5),
                         stepping = "bla")
    for stepping in ("", "bla", None):
        yield testStepping, stepping

def testParameterNumericalDisplayValues():
    dv = {1: 'one', 2: 'two', 3: 'three'}
    p = Parameter("testpar", 1, valueRange = (1, 5),
                     displayValues = dv)()
    for key, value in dv.iteritems():
        assert key in p.displayValues()
        assert p.displayValues(key) == value

def testParameterFloat():
    @raises(DecimalsError)
    def testDecimals(newDecimals):
        p = Parameter("testpar", 1.0, valueRange = (1, 5),
                         decimals = newDecimals)
    for value in ("", "bla", (1, 2), -1):
        yield testDecimals, value

def testParameterBaseCopy():
    p1 = Parameter(name = "p", value = "a", displayName = "displayname",
                   onValueUpdate = Dummy().dummyFunc)()
    p2 = p1.copy()
    assert isinstance(p1, ParameterBase)
    assert isinstance(p2, ParameterBase)
    assert p1 == p2
    p1.setValue("b")
    assert p1.value() != p2.value()
    p1.setDisplayName("q")
    assert p1.displayName() != p2.displayName()

def testParameterNumericalCopy():
    p1 = Parameter(name = "p", value = 1, displayName = "displayname",
                      valueRange = (1, 5), suffix = "suf", stepping = 1,
                      displayValues = {}, generator = NumberGenerator
                      )()
    p2 = p1.copy()
    assert isinstance(p1, ParameterNumerical)
    assert isinstance(p2, ParameterNumerical)
    assert p1 == p2
    p1.setValueRange((2, 3))
    assert p1.valueRange() != p2.valueRange()
    p1.setSuffix("suv")
    assert p1.suffix() != p2.suffix()
    p1.setStepping(2)
    assert p1.stepping() != p2.stepping()
    p1.setDisplayValues({1: "suv"})
    assert p1.displayValues() != p2.displayValues()
    p1.setGenerator(RandomUniform)
    assert p1.generator() != p2.generator()

def testParameterFloatCopy():
    p1 = Parameter(name = "p", value = 1.0, displayName = "displayname",
                      valueRange = (1, 5), suffix = "suf", stepping = 1,
                      displayValues = {}, generator = NumberGenerator,
                      decimals = 2
                      )()
    p2 = p1.copy()
    assert isinstance(p1, ParameterBase)
    assert isinstance(p2, ParameterBase)
    assert p1 == p2
    p1.setDecimals(4)
    assert p1.decimals() != p2.decimals()

def testParameterCompare():
    def compareParameters(pType, kwargs):
        p1 = Parameter(**kwargs)()
        p2 = Parameter(**kwargs)()
        assert isinstance(p1, pType)
        assert p1 == p2
    for pType, kwargs in (
            (ParameterBase,
                dict(name = "p", value = "1")),
            (ParameterBase,
                dict(name = "p", value = "a", displayName = "displayname")),
            (ParameterBase,
                dict(name = "p", value = "12", displayName = "displayname",
                     description = "description")),
            (ParameterNumerical,
                dict(name = "p", value = 1, valueRange = (1, 5),
                     onValueUpdate = Dummy().dummyFunc)),
            (ParameterNumerical,
                dict(name = "p", value = 1.0, valueRange = (1, 5),
                     suffix = "suf")),
            (ParameterNumerical,
                dict(name = "p", value = 1.0, valueRange = (1, 5),
                     suffix = "suf", stepping = 1.0)),
            (ParameterFloat,
                dict(name = "p", value = 1.0, valueRange = (1, 5))),
            (ParameterFloat,
                dict(name = "p", value = 1.0, valueRange = (1, 5),
                     decimals = 4)),
            (ParameterLog,
                dict(name = "p", value = 1.0, valueRange = (2., 4),
                     cls = ParameterLog)),
            ):
        yield compareParameters, pType, kwargs

def testParameterSerialize():
    import pickle
    import utils.pickleinstancemethods
    pType = factory(name = "radius", value = 3.4, valueRange = (1,5), decimals = 3)
    p = pType()
    p.setOnValueUpdate(Dummy().dummyFunc)
    data = pickle.dumps(p)
    p2 = pickle.loads(data)
    assert p == p2

# vim: set ts=4 sts=4 sw=4 tw=0:
