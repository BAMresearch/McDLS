# -*- coding: utf-8 -*-
# parameter.py

from math import log10 as math_log10
from math import fabs as math_fabs
from utils import isString, isNumber, isList
from inspect import getmembers

class ParameterError(StandardError):
    pass

class Parameter(object):
    """Base class for algorithm parameters."""
    _name = None
    _displayName = None
    _defaultValue = None

    def __init__(self, name, defaultValue,
                 displayName = None, description = None):
        assert defaultValue is not None, \
                "Filter parameter default value is mandatory!"
        self._defaultValue = defaultValue
        assert isString(name) and len(name) > 0, \
                "Filter parameter name is mandatory!"
        assert name.find(" ") < 0, \
                "Parameter name must not contain white space!"
        self._name = name
        if (not isString(displayName) or
            len(displayName) <= 0):
            displayName = self.name
        self._displayName = displayName
        if isString(description):
            self.__doc__ = description

    @property
    def name(self):
        return self._name

    @property
    def displayName(self):
        return self._displayName

    @property
    def defaultValue(self):
        return self._defaultValue

    def __str__(self):
        return str((self.name, self.displayName, self.defaultValue))

    def __eq__(self, other):
        for name, val in getmembers(self):
            if name.startswith("_"):
                continue
            if not hasattr(other, name):
                return False
            if getattr(self, name) != getattr(other, name):
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

class ParameterNumerical(Parameter):
    _valueRange = None
    _suffix = None
    _stepping = None
    _displayValues = None

    def __init__(self, name, defaultValue, valueRange,
                 suffix = None, stepping = None,
                 displayValues = None, **kwargs):
        Parameter.__init__(self, name, defaultValue, **kwargs)
        assert isNumber(defaultValue), \
                "The default value has to be numerical!"
        assert valueRange is not None, \
                "The value range is mandatory for a numerical parameter!"
        assert len(valueRange) == 2
        assert isNumber(valueRange[0])
        assert isNumber(valueRange[1])
        self._valueRange = tuple(valueRange)
        if isString(suffix) and len(suffix) > 0:
            self._suffix = suffix
        if isNumber(stepping):
            self._stepping = stepping
        self._displayValues = dict()
        if isList(displayValues) and all([isList(x) for x in displayValues]):
            for valueList in displayValues:
                self._displayValues[valueList[0]] = valueList[1]
                self._displayValues[valueList[1]] = valueList[0]

    @property
    def valueRange(self):
        return self._valueRange

    @property
    def suffix(self):
        return self._suffix

    @property
    def stepping(self):
        return self._stepping

    @property
    def displayValues(self):
        return self._displayValues

    def __str__(self):
        return Parameter.__str__(self) + str((
                self.valueRange, self.suffix, self.stepping))

class ParameterFloat(ParameterNumerical):
    _decimals = None

    def __init__(self, name, defaultValue, valueRange, decimals = None,
                 **kwargs):
        ParameterNumerical.__init__(self, name, defaultValue, valueRange,
                                    **kwargs)
        if isNumber(decimals):
            self._decimals = int(decimals)
        else:
            start, end = self.valueRange
            self._decimals = int(round(math_log10(math_fabs(end - start))))
        if self._decimals <= 0:
            self._decimals = 1

    @property
    def decimals(self):
        return self._decimals

    def __str__(self):
        return ParameterNumerical.__str__(self) + str((self.decimals,))

class ParameterLog(ParameterFloat):
    pass

# vim: set ts=4 sts=4 sw=4 tw=0:
