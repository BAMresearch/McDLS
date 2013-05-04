# -*- coding: utf-8 -*-
# parameter_test.py

from parameter import (Parameter, ParameterNumerical, ParameterFloat,
    ParameterNameError, DefaultValueError, ValueRangeError, SuffixError,
    SteppingError, DecimalsError, DisplayValuesError)
from nose.tools import raises

def testParameterName():
    @raises(ParameterNameError)
    def testName(newName):
        class TestPar(Parameter):
            name = newName
            defaultValue = 0
        p = TestPar()
    for name in (None, "", 1.3, 0):
        yield testName, name

@raises(DefaultValueError)
def testParameterDefaultValue():
    class TestPar(Parameter):
        defaultValue = None
    p = TestPar()

def testParameterNumerical():
    class TestPar(ParameterNumerical):
        name = "testpar"
        defaultValue = 3
        valueRange = (1, 5)
        suffix = "mm"
        stepping = 1
    p = TestPar()
    assert p.defaultValue == 3
    assert p.valueRange == (1, 5)
    assert p.suffix == "mm"
    assert p.stepping == 1
    assert p.displayValues is None

def testParameterNumericalValueRange():
    @raises(ValueRangeError)
    def testValueRange(newRange):
        class TestPar(ParameterNumerical):
            name = "testpar"
            defaultValue = 1.0
            valueRange = newRange
        p = TestPar()
    for valueRange in (None, (None, 1), (1, None), (None, None), (1, 2, 3),
                       "", ("", ), (1, ""), ("", 1), ("", 1.0), (1.0, "")):
        yield testValueRange, valueRange

def testParameterNumericalSuffix():
    @raises(SuffixError)
    def testSuffix(newSuffix):
        class TestPar(ParameterNumerical):
            name = "testpar"
            defaultValue = 1.0
            valueRange = (1, 5)
            suffix = newSuffix
        p = TestPar()
    for suffix in ("", 1, 1.0):
        yield testSuffix, suffix

def testParameterNumericalStepping():
    @raises(SteppingError)
    def testStepping(stepping):
        class TestPar(ParameterNumerical):
            name = "testpar"
            defaultValue = 1.0
            valueRange = (1, 5)
            stepping = "bla"
        p = TestPar()
    for stepping in ("", "bla", None):
        yield testStepping, stepping

def testParameterNumericalDisplayValues():
    dv = {1: 'one', 2: 'two', 3: 'three'}
    class TestPar(ParameterNumerical):
        name = "testpar"
        defaultValue = 1.0
        valueRange = (1, 5)
        displayValues = dv
    p = TestPar()
    for key, value in dv.iteritems():
        assert key in p.displayValues
        assert p.displayValues[key] == value

def testParameterFloat():
    @raises(DecimalsError)
    def testDecimals(newDecimals):
        class TestPar(ParameterFloat):
            name = "testpar"
            defaultValue = 1.0
            valueRange = (1, 5)
            decimals = newDecimals
        p = TestPar()
    for value in ("", "bla", (1, 2), 0, -1):
        yield testDecimals, value

# vim: set ts=4 sts=4 sw=4 tw=0:
