# -*- coding: utf-8 -*-
# parameter.py

from math import log10 as math_log10
from math import fabs as math_fabs
from utils import isString, isNumber, isList, isMap
from inspect import getmembers

class ParameterError(StandardError):
    pass

class DefaultValueError(ParameterError):
    pass

class ParameterNameError(ParameterError):
    pass

class ValueRangeError(ParameterError):
    pass

class SuffixError(ParameterError):
    pass

class SteppingError(ParameterError):
    pass

class DecimalsError(ParameterError):
    pass

class DisplayValuesError(ParameterError):
    pass

def testfor(condition, exception, errorMessage):
    if not __debug__:
        return
    if not condition:
        raise exception(errorMessage)

class Parameter(object):
    """Base class for algorithm parameters providing additional
    information to ease automated GUI building."""
    # class variables, all objects of the same type shared these values
    name = None
    displayName = None
    defaultValue = None
    # instance variables, every object has its own value
    _value = None

    @classmethod
    def __str__(cls):
        return str((cls.name, cls.displayName, cls.defaultValue))

    def __init__(self, value = None):
        # Checking attributes which have to be set in a class of this type.
        testfor(self.defaultValue is not None,
                DefaultValueError, "Default value is mandatory!")
        testfor(isString(self.name) and len(self.name) > 0,
                ParameterNameError, "A name is mandatory!")
        testfor(self.name.find(" ") < 0,
                ParameterNameError, "A name must not contain white space!")
        if (not isString(self.displayName) or len(self.displayName) <= 0):
            self.displayName = self.name
        self.value = self.defaultValue
        if isNumber(value):
            self.value = value

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, newValue):
        self._value = newValue

class ParameterNumerical(Parameter):
    # class variables, see Parameter
    valueRange = None
    suffix = None
    stepping = None
    displayValues = None # dict maps values to text being displayed instead
    # instance variables, see Parameter
    _error = None

    @staticmethod
    def assertRange(newRange):
        testfor(isList(newRange), ValueRangeError,
                "A value range is mandatory for a numerical parameter!")
        testfor(len(newRange) == 2, ValueRangeError,
                "A value range has to consist of two values!")
        testfor(all([isNumber(v) for v in newRange]), ValueRangeError,
                "A value range has to consist of numbers only!")

    @classmethod
    def __str__(cls):
        return Parameter.__str__() + str((
                cls.valueRange, cls.suffix, cls.stepping))

    def __init__(self, *args, **kwargs):
        Parameter.__init__(self, *args, **kwargs)
        testfor(isNumber(self.defaultValue), DefaultValueError,
                "A default value has to be numerical!")
        self.assertRange(self.valueRange)
        if self.suffix is not None:
            testfor(isString(self.suffix) and len(self.suffix) > 0,
                    SuffixError, "Parameter suffix has to be some text!")
        if self.stepping is not None:
            testfor(isNumber(self.stepping),
                    SteppingError, "Parameter has to be a number!")
        if self.displayValues is not None:
            testfor(isMap(self.displayValues), DisplayValuesError,
                    "Expected a display value mapping of numbers to text!")
            testfor(all([isNumber(v) for v in self.displayValues.iterkeys()]),
                DisplayValuesError, "Display value keys have to be numbers!")
            testfor(all([isString(s) for s in self.displayValues.itervalues()]),
                DisplayValuesError, "Display values have to be text!")

    @property
    def valueRange(self):
        return self.valueRange

    @valueRange.setter
    def valueRange(self, newRange):
        cls.assertRange(newRange)
        self.valueRange = newRange

    @property
    def error(self):
        return self._error

    @error.setter
    def error(self, newError):
        msg = "A value error has to be of the same data type than the value!"
        if isList(newError):
            testfor(all([isinstance(e, type(self.value))
                        for e in newError]), ParameterErrorError, msg)
        else:
            testfor(isinstance(newError, type(self.value)),
                    ParameterErrorError, msg)
        self._error = newError

class ParameterFloat(ParameterNumerical):
    decimals = None

    @classmethod
    def __str__(cls):
        return ParameterNumerical.__str__() + str((cls.decimals,))

    def __init__(self, *args, **kwargs):
        ParameterNumerical.__init__(self, *args, **kwargs)
        if self.decimals is not None:
            testfor(isNumber(self.decimals) and self.decimals > 0,
                    DecimalsError, "Parameter decimals has to be a number!")
            self.decimals = int(self.decimals)
        else:
            start, end = self.valueRange
            self.decimals = int(round(math_log10(math_fabs(end - start))))
        if self.decimals <= 0:
            self.decimals = 1

class ParameterLog(ParameterFloat):
    pass

# vim: set ts=4 sts=4 sw=4 tw=0:
