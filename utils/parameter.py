# -*- coding: utf-8 -*-
# parameter.py

from math import log10 as math_log10
from math import fabs as math_fabs
from inspect import getmembers
import logging
import numpy
import sys
from utils import isString, isNumber, isList, isMap
from utils.mixedmethod import mixedmethod
from utils.classproperty import classproperty
from numbergenerator import NumberGenerator, RandomUniform

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

class ParameterGeneratorError(ParameterError):
    pass

def testfor(condition, exception, errorMessage):
    if __debug__ and not condition:
        raise exception(errorMessage)

# we have to differentiate class variables and instance variables
# class variables: names, suffix, stepping, datatype
# -> valueRange? general definition range in type
# -> user provided subrange for actual calculations in instance
# only the instance has values and ranges
# -> everything which differs for other scattering contributions

class Parameter(object):
    """Base class for algorithm parameters providing additional
    information to ease automated GUI building."""

    _name = None
    _displayName = None
    _value = None

    @staticmethod
    def setup(cls, name, value, displayName = None):
        def assertName(newName):
            testfor(isString(newName) and len(newName) > 0,
                    ParameterNameError, "A name is mandatory!")
            testfor(newName.find(" ") < 0,
                    ParameterNameError, "A name must not contain white space!")
        assertName(name)
        newParam = type(name.title()+"Parameter", (cls,), dict())
        newParam._name = name
        newParam.setValue(value)
        newParam.setDisplayName(displayName)
        return newParam

    @classmethod
    def make(cls, *args, **kwargs):
        return cls.setup(cls, *args, **kwargs)

    @mixedmethod
    def setValue(selforcls, newValue):
        testfor(newValue is not None,
                DefaultValueError, "Default value is mandatory!")
        selforcls._value = newValue

    @mixedmethod
    def setDisplayName(selforcls, newName):
        if (not isString(newName) or len(newName) <= 0):
            newName = selforcls.name()
        selforcls._displayName = newName

    @mixedmethod
    def name(self):
        return self._name

    @mixedmethod
    def value(self):
        return self._value

    @mixedmethod
    def displayName(self):
        return self._displayName

    @classproperty
    @classmethod
    def dtype(cls):
        return str

    def __str__(self):
        return "{0}: {1}".format(
                self.displayName(), self.value())

class ParameterNumerical(Parameter):
    _valueRange = None
    _suffix = None
    _stepping = None
    _displayValues = None # dict maps values to text being displayed instead
    _generator = None

    @staticmethod
    def setup(cls, name, value, displayName = None, valueRange = None,
             suffix = None, stepping = None, displayValues = None,
             generator = None):
        newParam = Parameter.setup(cls, name, value, displayName)
        newParam.setValueRange(valueRange)
        newParam.setSuffix(suffix)
        newParam.setStepping(stepping)
        newParam.setDisplayValues(displayValues)
        newParam.setGenerator(generator)
        return newParam

    @mixedmethod
    def setValue(selforcls, newValue):
        testfor(isNumber(newValue), DefaultValueError,
                "A default value has to be numerical!")
        super(ParameterNumerical, selforcls).setValue(newValue)

    @mixedmethod
    def setValueRange(selforcls, newRange):
        testfor(isList(newRange), ValueRangeError,
                "A value range is mandatory for a numerical parameter!")
        testfor(len(newRange) == 2, ValueRangeError,
                "A value range has to consist of two values!")
        testfor(all([isNumber(v) for v in newRange]), ValueRangeError,
                "A value range has to consist of numbers only!")
        selforcls._valueRange = minVal, maxVal = min(newRange), max(newRange)
        if selforcls._value < minVal:
            selforcls._value = minVal
        if selforcls._value > maxVal:
            selforcls._value = maxVal

    @mixedmethod
    def setSuffix(selforcls, newSuffix):
        if newSuffix is None:
            return
        testfor(isString(newSuffix) and len(newSuffix) > 0,
                SuffixError, "Parameter suffix has to be some text!")
        selforcls._suffix = newSuffix

    @mixedmethod
    def setStepping(selforcls, newStepping):
        if newStepping is None:
            return
        testfor(isNumber(newStepping),
                SteppingError, "Parameter has to be a number!")
        selforcls._stepping = newStepping

    @mixedmethod
    def setDisplayValues(selforcls, newDisplayValues):
        if newDisplayValues is None:
            return
        testfor(isMap(newDisplayValues), DisplayValuesError,
                "Expected a display value mapping of numbers to text!")
        testfor(all([isNumber(v) for v in newDisplayValues.iterkeys()]),
            DisplayValuesError, "Display value keys have to be numbers!")
        testfor(all([isString(s) for s in newDisplayValues.itervalues()]),
            DisplayValuesError, "Display values have to be text!")
        selforcls._displayValues = newDisplayValues

    @mixedmethod
    def setGenerator(selforcls, newGenerator):
        if isinstance(newGenerator, type):
            testfor(issubclass(newGenerator, NumberGenerator),
                    ParameterGeneratorError, "NumberGenerator type expected!")
        else:
            newGenerator = RandomUniform
        selforcls._generator = newGenerator
        logging.info("Parameter {0} uses {1} distribution."
                     .format(selforcls._name, newGenerator.__name__))

    @mixedmethod
    def valueRange(self, index = None):
        if not isNumber(index) or index not in range(len(self._valueRange)):
            return self._valueRange
        else:
            return self._valueRange[index]

    @mixedmethod
    def suffix(self):
        return self._suffix

    @mixedmethod
    def stepping(self):
        return self._stepping

    @mixedmethod
    def displayValues(self, key = None):
        if key is None:
            return self._displayValues
        else:
            return self._displayValues.get(key, None)

    @mixedmethod
    def generator(self):
        return self._generator

    @classproperty
    @classmethod
    def dtype(cls):
        return int

    def __str__(self):
        return (Parameter.__str__(self) + " in [{0}, {1}]{2}, {3} steps"
                .format(*(self.valueRange() + (self.suffix(),
                                               self.stepping()))))

    def generate(self, lower = None, upper = None, count = 1):
        raise NotImplementedError

class ParameterFloat(ParameterNumerical):
    _decimals = None

    @staticmethod
    def setup(cls, name, value, displayName = None, valueRange = None,
             suffix = None, stepping = None, displayValues = None,
             generator = None, decimals = None):
        newParam = ParameterNumerical.setup(cls, name, value, displayName,
                    valueRange, suffix, stepping, displayValues, generator)
        newParam = type(newParam._name, (newParam,), dict())
        newParam.setDecimals(decimals)
        return newParam

    @mixedmethod
    def setDecimals(selforcls, newDecimals):
        if newDecimals is not None:
            testfor(isNumber(newDecimals) and newDecimals > 0,
                    DecimalsError, "Parameter decimals has to be a number!")
        else:
            start, end = selforcls._valueRange
            newDecimals = round(math_log10(math_fabs(end - start)))
        newDecimals = max(newDecimals, 1)
        newDecimals = min(newDecimals, sys.float_info.max_10_exp)
        selforcls._decimals = int(newDecimals)

    @mixedmethod
    def decimals(self):
        return self._decimals

    @classproperty
    @classmethod
    def dtype(cls):
        return float

    def __str__(self):
        return (ParameterNumerical.__str__(self) +
                ", {0} decimals".format(self.decimals()))

    def generate(self, lower = None, upper = None, count = 1):
        """Returns a list of valid parameter values within given bounds.
        Accepts vectors of individual bounds for lower and upper limit.
        This allows for inequality parameter constraints.
        """
        # works with vectors of multiple bounds too
        vRange = self.valueRange()
        if lower is None:
            lower = vRange[0]
        if upper is None:
            upper = vRange[1]
        vRange = (numpy.maximum(vRange[0], lower),
                  numpy.minimum(vRange[1], upper))
        if isList(vRange[0]) and isList(vRange[1]):
            assert len(vRange[0]) == len(vRange[1]), \
                "Provided value range is unsymmetrical!"
        try: # update count to length of provided bound vectors
            count = max(count, min([len(x) for x in vRange]))
        except:
            pass
        values = self.generator().get(count)
        # scale numbers to requested range
        return values * (vRange[1] - vRange[0]) + vRange[0]

class ParameterLog(ParameterFloat):
    pass

# vim: set ts=4 sts=4 sw=4 tw=0:
