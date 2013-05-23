# -*- coding: utf-8 -*-
# parameter_test.py

from parameter import (Parameter, ParameterNumerical, ParameterFloat,
    ParameterNameError, DefaultValueError, ValueRangeError, SuffixError,
    SteppingError, DecimalsError, DisplayValuesError)
from nose.tools import raises

def testParameterName():
    @raises(ParameterNameError)
    def testName(newName):
        p = Parameter.make(name, 0)
    for name in (None, "", 1.3, 0):
        yield testName, name

@raises(TypeError)
def testParameterDefaultValue1():
    p = Parameter.make("testpar")

@raises(DefaultValueError)
def testParameterDefaultValue2():
    p = Parameter.make("testpar", None)

def testParameterNumerical():
    ptype = ParameterNumerical.make("testpar", 3, valueRange = (1, 5),
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
        p = ParameterNumerical.make("testpar", 1.0, valueRange = newRange)
    for valueRange in (None, (None, 1), (1, None), (None, None), (1, 2, 3),
                       "", ("", ), (1, ""), ("", 1), ("", 1.0), (1.0, "")):
        yield testValueRange, valueRange

def testParameterNumericalSuffix():
    @raises(SuffixError)
    def testSuffix(newSuffix):
        p = ParameterNumerical.make("testpar", 1.0, valueRange = (1, 5),
                                    suffix = newSuffix)
    for suffix in ("", 1, 1.0):
        yield testSuffix, suffix

def testParameterNumericalStepping():
    @raises(SteppingError)
    def testStepping(stepping):
        p = ParameterNumerical.make("testpar", 1.0, valueRange = (1, 5),
                                    stepping = "bla")
    for stepping in ("", "bla", None):
        yield testStepping, stepping

def testParameterNumericalDisplayValues():
    dv = {1: 'one', 2: 'two', 3: 'three'}
    p = ParameterNumerical.make("testpar", 1.0, valueRange = (1, 5),
                                displayValues = dv)()
    for key, value in dv.iteritems():
        assert key in p.displayValues()
        assert p.displayValues(key) == value

def testParameterFloat():
    @raises(DecimalsError)
    def testDecimals(newDecimals):
        p = ParameterFloat.make("testpar", 1.0, valueRange = (1, 5),
                                decimals = newDecimals)
    for value in ("", "bla", (1, 2), 0, -1):
        yield testDecimals, value

# vim: set ts=4 sts=4 sw=4 tw=0:
