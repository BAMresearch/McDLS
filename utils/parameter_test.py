# -*- coding: utf-8 -*-
# parameter_test.py

from parameter import Parameter, ParameterNumerical, ParameterFloat
from nose.tools import raises

def testParameter():
    p = Parameter("testpar", "bla")
    assert p.displayName == "testpar"
    assert p.defaultValue == "bla"

    p = Parameter("testpar", "bla", "Test Par", "test parameter")
    assert p.defaultValue == "bla"
    assert p.displayName == "Test Par"
    assert p.__doc__ == "test parameter"

def testParameterName():
    @raises(AssertionError)
    def testName(name):
        p = Parameter(name, "")
    for name in (None, "", 1.3, 0):
        yield testName, name

@raises(AssertionError)
def testParameterDefaultValue():
    p = Parameter("testpar", None)

def testParameterNumerical():
    p = ParameterNumerical("testpar", 3, (1, 5), "mm", 1)
    assert p.defaultValue == 3
    assert p.valueRange == (1, 5)
    assert p.suffix == "mm"
    assert p.stepping == 1
    assert p.displayValues == dict()

#@raises(AssertionError)
#def testParameterNumericalDefaultValue():
#    p = ParameterNumerical("testpar", "bla", (1, 5))

def testParameterNumericalValueRange():
    @raises(AssertionError)
    def testValueRange(valueRange):
        p = ParameterNumerical("testpar", 1.0, valueRange)
    for valueRange in (None, (None, 1), (1, None), (None, None), (1, 2, 3),
                       "", ("", ), (1, ""), ("", 1), ("", 1.0), (1.0, "")):
        yield testValueRange, valueRange

def testParameterNumericalSuffix():
    def testSuffix(suffix):
        p = ParameterNumerical("testpar", 1.0, (1, 5), suffix)
        assert p.suffix == None
    for suffix in ("", None, 1, 1.0):
        yield testSuffix, suffix

def testParameterNumericalStepping():
    def testStepping(stepping):
        p = ParameterNumerical("testpar", 1.0, (1, 5), stepping = "bla")
        assert p.stepping == None
    for stepping in ("", "bla", None):
        yield testStepping, stepping

def testParameterNumericalDisplayValues():
    dv = ((1, 'one'), (2, 'two'), (3, 'three'))
    p = ParameterNumerical("testpar", 1.0, (1, 5),
            displayValues = dv)
    for v in dv:
        assert v[0] in p.displayValues
        assert v[1] in p.displayValues

def testParameterFloat():
    def testDecimals(value):
        p = ParameterFloat("testpar", 1.0, (1, 5), decimals = value)
        assert p.decimals == 1 # derived from value range
    for value in (None, "", "bla", (1, 2), 0, -1, 1.5):
        yield testDecimals, value

def testParameterCompare():
    for pType, args in ((Parameter, ["p", 1.0]),
                        (Parameter, ["p", 1.0, "displayname"]),
                        (Parameter, ["p", 1.0, "displayname", "description"]),
                        (ParameterNumerical, ["p", 1.0, (1, 5)]),
                        (ParameterNumerical, ["p", 1.0, (1, 5), "suf"]),
                        (ParameterNumerical, ["p", 1.0, (1, 5), "suf", 1.0]),
                        (ParameterFloat, ["p", 1.0, (1, 5)]),
                        (ParameterFloat, ["p", 1.0, (1, 5), 4]),
                        ):
        p1 = pType(*args)
        p2 = pType(*args)
        assert p1 == p2

# vim: set ts=4 sts=4 sw=4 tw=0:
