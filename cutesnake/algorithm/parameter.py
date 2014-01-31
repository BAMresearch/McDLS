# -*- coding: utf-8 -*-
# parameter.py

"""
This module defines a generic parameter class for algorithms.
It contains meta information which allows for automated UI building.
Create sub classes by calling factory() in this module.
It creates a new sub class type which inherits ParameterBase::

>>> from parameter import factory as paramFactory
>>> ParamType = paramFactory("radius", 1.3, valueRange = (0, 2))

Created a new type RadiusParameter:

>>> print(ParamType)
<class 'parameter.RadiusParameter'>

Using methods on instances work as usual:

>>> p = ParamType()
>>> p.name()
'radius'
>>> p.value()
1.3

Update the instance:

>>> p.setValue(2.4)
>>> p.value()
2.4

Changing class default:
>>> ParamType.setValue(3.5)
>>> ParamType.value()
3.5

Existing instance keep their values:
>>> p.value()
2.4

New instances get the updated defaults:
>>> q = ParamType()
>>> q.value()
3.5

Parameter attributes are accessible on type/class as well as on the
instance. Updating an attribute of an instance changes just that
individual instance whereas updating an attribute of the type changes
that attribute in general for all new instances to be created
which is behaves like a default value.
"""

from math import log10 as math_log10
from math import fabs as math_fabs
from inspect import getmembers
import numpy
import sys
import logging
from cutesnake.utils import isString, isNumber, isList, isMap
from cutesnake.utils.tests import testfor, assertName
from cutesnake.utils.mixedmethod import mixedmethod
from cutesnake.utils.classproperty import classproperty
from numbergenerator import NumberGenerator, RandomUniform

def generateValues(numberGenerator, defaultRange, lower, upper, count):
    # works with vectors of multiple bounds too
    vRange = defaultRange
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
    values = numberGenerator.get(count)
    # scale numbers to requested range
    return values * (vRange[1] - vRange[0]) + vRange[0]

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

def _makeGetter(varName):
    def getter(selforcls):
        return getattr(selforcls, varName)
    return mixedmethod(getter)

def _makeSetter(varName):
    def setter(selforcls, value):
        setattr(selforcls, varName, value)
    return mixedmethod(setter)

def _setterName(attrName):
    return "set" + attrName[0].upper() + attrName[1:]

class ParameterBase(object):
    """Base class for algorithm parameters providing additional
    information to ease automated GUI building."""

    # Be able to manage attributes programmatically also for derived classes
    # while maintaining the order of attributes which is not preserved by
    # pythons __dict__. Initialization of *decimals* has to occur after
    # *valueRange*. Not sure, if this is a good idea but it reduces
    # repetition/code drastically.
    @classmethod
    def setAttributes(cls, dictionary, *names, **namesAndValues):
        """Sets an *ordered* list of attributes.
        Initializes the private variable to None and sets a default getter
        method for each name provided. Additionally, sets *attributeNames* to
        return all attribute names."""
        names += tuple(namesAndValues.keys())
        for attrName in names:
            varName = "_"+attrName
            # set the class variable
            dictionary[varName] = namesAndValues.get(attrName, None)
            # default getter for the class var
            dictionary[attrName] = _makeGetter(varName)
            if attrName != "name": # do not create setName() method
                dictionary[_setterName(attrName)] = _makeSetter(varName)
        try:
            # prepend attribute names of the base class
            names = cls.attributeNames() + names
        except:
            pass
        # sets the ordered names of attributes
        dictionary["_attributeNames"] = names

    setAttributes.__func__(None, locals(), "name", "value", "displayName")

    @classmethod
    def attributeNames(cls):
        """Returns an ordered list of attribute names considering multiple
        inheritance and maintaining its order."""
        mergedAttrNames = []
        for baseCls in reversed(cls.__mro__):
            if not hasattr(baseCls, "_attributeNames"):
                continue
            mergedAttrNames += [attrName
                                for attrName in baseCls._attributeNames
                                if attrName not in mergedAttrNames]
        return mergedAttrNames

    @mixedmethod
    def attributes(selforcls):
        """Returns a dictionary with <key, value> pairs of all attributes and
        their values in this type or instance.
        Helps to avoid having explicit long argument lists"""
        base = selforcls
        if isinstance(base, object):
            base = type(base)
        base = base.__mro__[1] # Parameter classes have only one direct subclass
        # store the direct base class for duplication later
        res = dict(cls = base, description = selforcls.__doc__)
        for name in selforcls.attributeNames():
            # if this throws an exception, there is a bug
            res[name] = getattr(selforcls, name)()
        return res

    @classmethod
    def factory(cls, **kwargs):
        """Returns this type with attribute values initialized in the ordering
        provided by attributeNames()"""
        for key, value in kwargs.iteritems():
            # set the attributes for which we find setters
            # the setter may raise exceptions for invalid data
            setter = getattr(cls, _setterName(key), None)
            if setter is not None: # key exists, check for method?
                setter(value)
        return cls

    @classmethod
    def copy(cls):
        attr = cls.attributes()
        other = attr.pop("cls") # remove duplicate first argument of factory()
        if not issubclass(other, cls):
            other = cls
        return other.factory(**attr)()

    @classmethod
    def get(self, key, default = None):
        """metagetter to get an attribute parameter"""
        if key in self.attributeNames():
            getterFunc = getattr(self, key, default)
            return getterFunc()
        else:
            logging.warning(
                    "parameter {n} attribute {k} not understood in get"
                    .format(n = self.name(), k = key))

    @classmethod
    def set(self, key, value):
        """metasetter to set an attribute value"""
        if key in self.attributeNames():
            setterFunc = getattr(self, _setterName(key))
            setterFunc(value)
        else:
            logging.warning(
                    "parameter {n} attribute {k} not found in set"
                    .format(n = self.name(), k = key))

    @classmethod
    def setName(cls, name):
        """Changing the name is allowed for the class/type only,
        not for instances."""
        assertName(name, ParameterNameError)
        cls._name = name

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

    @classproperty
    @classmethod
    def dtype(cls):
        return str

    @classmethod
    def isDataType(cls, value):
        return isinstance(value, cls.dtype)

    def __str__(self):
        return "{0}: {1}".format(
                self.displayName(), self.value())

    def __eq__(self, other):
        if self.dtype != other.dtype:
            return False
        try:
            return self.attributes() == other.attributes()
        except:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def generate(self, lower = None, upper = None, count = 1):
        """Returns a list of valid parameter values within given bounds.
        Accepts vectors of individual bounds for lower and upper limit.
        This allows for inequality parameter constraints.
        lower, upper: arrays for lower and upper bounds
        """
        raise NotImplementedError

class ParameterBoolean(ParameterBase):
    @classproperty
    @classmethod
    def dtype(cls):
        return bool

class ParameterString(ParameterBase):
    """
    String-based parameter class. The default value should be the first
    item in the _valueRange list.
    """
    ParameterBase.setAttributes(locals(), "valueRange")

    @classmethod
    def setValueRange(self, newRange):
        testfor(isList(newRange), ValueRangeError,
                "A value range for a string type parameter has to be a list!")
        testfor(all([isString(v) for v in newRange]), ValueRangeError,
                "A value range for a string has to be a list of strings!")
        self._valueRange = newRange
        if not (self.value() in self.valueRange()):
            # where are the default values?
            self.setValue(self.valueRange()[0])

    @mixedmethod
    def valueRange(self):
        if self._valueRange is None:
            return ()
        return self._valueRange

    @classproperty
    @classmethod
    def dtype(cls):
        return str

    @classmethod
    def isDataType(cls, value):
        return isString(value) 

class ParameterNumerical(ParameterBase):
    # defines attributes for this parameter type and creates default
    # getter/setter for them. For specialized versions they can be
    # overridden as usual.
    ParameterBase.setAttributes(locals(), "valueRange", "suffix",
                  "stepping", "displayValues", "generator")

    @mixedmethod
    def setValue(selforcls, newValue):
        testfor(isNumber(newValue), DefaultValueError,
                "A value has to be numerical!")
        super(ParameterNumerical, selforcls).setValue(newValue)

    @mixedmethod
    def setValueRange(selforcls, newRange):
        testfor(isList(newRange), ValueRangeError,
                "A value range is mandatory for a numerical parameter!")
        testfor(len(newRange) == 2, ValueRangeError,
                "A value range has to consist of two values!")
        testfor(all([isNumber(v) for v in newRange]), ValueRangeError,
                "A value range has to consist of numbers only!")
        minVal, maxVal = min(newRange), max(newRange)
        #minVal = max(minVal, -sys.float_info.max)
        #maxVal = min(maxVal,  sys.float_info.max)
        minVal = max(minVal, -1e200)
        maxVal = min(maxVal,  1e200)
        selforcls._valueRange = minVal, maxVal
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
        # TODO: also add reverse lookup
        selforcls._displayValues = newDisplayValues

    @mixedmethod
    def setGenerator(selforcls, newGenerator):
        if isinstance(newGenerator, type):
            testfor(issubclass(newGenerator, NumberGenerator),
                    ParameterGeneratorError, "NumberGenerator type expected!")
        else:
            newGenerator = RandomUniform
        selforcls._generator = newGenerator
#        logging.info("Parameter {0} uses {1} distribution."
#                     .format(selforcls._name, newGenerator.__name__))

    @mixedmethod
    def valueRange(selforcls):
        if selforcls._valueRange is None:
            return (None, None)
        return selforcls._valueRange

    @mixedmethod
    def min(selforcls):
        return selforcls.valueRange()[0]

    @mixedmethod
    def max(selforcls):
        return selforcls.valueRange()[1]

    @mixedmethod
    def displayValues(selforcls, key = None, default = None):
        if key is None:
            return selforcls._displayValues
        else:
            return selforcls._displayValues.get(key, default)

    @classproperty
    @classmethod
    def dtype(cls):
        return int

    @classmethod
    def isDataType(cls, value):
        """ParameterNumerical is a fallback for all number not being float."""
        return isNumber(value) and not isinstance(value, float)

    def __str__(self):
        return (ParameterBase.__str__(self) + " in [{0}, {1}]{2}, {3} steps"
                .format(*(self.valueRange() + (self.suffix(),
                                               self.stepping()))))

    def generate(self, lower = None, upper = None, count = 1):
        return generateValues(self.generator(), self.valueRange(),
                                                lower, upper, count
                                ).astype(self.dtype)

class ParameterFloat(ParameterNumerical):
    ParameterNumerical.setAttributes(locals(), "decimals")

    @mixedmethod
    def setDecimals(selforcls, newDecimals):
        if newDecimals is not None:
            testfor(isNumber(newDecimals) and newDecimals >= 0, DecimalsError,
                    "Parameter decimals has to be a positive number!")
        else:
            start, end = selforcls._valueRange
            newDecimals = round(math_log10(math_fabs(end - start)))
        newDecimals = max(newDecimals, 0)
        newDecimals = min(newDecimals, sys.float_info.max_10_exp)
        selforcls._decimals = int(newDecimals)

    @classproperty
    @classmethod
    def dtype(cls):
        return float

    @classmethod
    def isDataType(cls, value):
        return isinstance(value, cls.dtype)

    def __str__(self):
        return (ParameterNumerical.__str__(self) +
                ", {0} decimals".format(self.decimals()))

    def generate(self, lower = None, upper = None, count = 1):
        return generateValues(self.generator(), self.valueRange(),
                                                lower, upper, count)

class ParameterLog(ParameterFloat):
    """Used to select an UI input widget with logarithmic behaviour."""
    pass

def factory(name, value, paramTypes = None, **kwargs):
    """
    Generates a new Parameter type derived from one of the predefined
    base classes choosen by the supplied value: Providing a string value
    results in a type derived from ParameterBase, providing an integer
    value produces a ParameterNumerical type and a float value results
    in a ParameterFloat type.
    Alternatively, a class type cls can be provided which is used as base
    class for the resulting Parameter class type. Make sure in this case,
    all attributes mandatory for this base type are provided too.

    *name*: short name of the new parameter without spaces
    *value*: default value from which the type is derived if cls is not given

    Optional arguments:
    *paramTypes*:  tuple of available parameter types instead of the default
    *cls*:         forces a certain Parameter type.
    *description*: Updates the __doc__ attribute. May be displayed in the UI
                   somewhere.
    """
    kwargs.update(name = name, value = value)
    name = kwargs.get("name", None)
    assertName(name, ParameterNameError)
    value = kwargs.get("value", None)
    cls = kwargs.pop("cls", None) # remove 'cls' keyword before forwarding
    if paramTypes is None:
        paramTypes = (ParameterBoolean, ParameterFloat,
                      ParameterNumerical, ParameterBase)
    if not (cls is not None and (
                (isinstance(cls, super) and
                 issubclass(cls.__thisclass__, ParameterBase)) or
                 issubclass(cls, ParameterBase))):
        for cls in paramTypes[:-1]:
            if cls.isDataType(value):
                break
        else:
            cls = paramTypes[-1] # ParameterBase usually
    # embed description as class documentation
    clsdict = dict()
    description = kwargs.get("description", None)
    if isString(description) and len(description) > 0:
        clsdict['__doc__'] = description
    # create a new class/type with given name and base class
    # translate works different for unicode strings:
    typeName = str(name.title()).translate(None, ' \t\n\r') + "Parameter"
    NewType = type(typeName, (cls,), clsdict)
    # set up the new class before return
    return NewType.factory(**kwargs)
    # creating new types here is a problem for pickle
    # it is done to be able to change/add class (type) attributes
    # without modifying the base classes
    # ATM: __doc__ + all class attribs get mod by Param.factory()

if __name__ == "__main__":
    import doctest
    doctest.testmod()

# vim: set ts=4 sts=4 sw=4 tw=0:
