# -*- coding: utf-8 -*-
# utils/units.py

u"""
Defines methods for using and manipulating units of variables. 
Some default magnitude-name dictionaries are provided, but the user can 
supply their own dictionary if required. Default unit to translate to 
must be set.
Required keyword arguments:

    - *magnitudedict*: a dictionary of magnitude - name pairs. Names must be
                       unicode strings.
    - *simagnitudename*: the si magnitude name.

Example usage: 

>>> rUnit = Length("nm")
>>> rUnit.siMagnitudeName
u'm'
>>> rUnit.displayMagnitudeName
u'nm'
>>> rUnit.magnitudeConversion
1e-09

or:
>>> rUnit.toSi(32)
3.2e-08

Selecting a default: 
>>> qUnit = ScatteringVector(u"cm⁻¹")
>>> qUnit.magnitudeConversion
100.0

"""

import logging
import collections
import numpy as np # For arrays
from numpy import pi
from utils.classproperty import classproperty

class Unit(object):
    _magnitudeMap = None
    _siMagnitudeName = None
    _displayMagnitudeName = None

    def __init__(self, magnitudeName = None):
        """Set up a unit of measurement. The provided magnitude name is used
        by default in the UI. Using SI units if unspecified."""
        if magnitudeName is None:
            magnitudeName = self.siMagnitudeName
        elif magnitudeName not in self.magnitudeMapping:
            logging.warning(u"Provided default magnitude name '{mn}' is not "
                            u"available: '{map}'!".format(
                            mn = magnitudeName, map = self.magnitudeMapping))
        self._displayMagnitudeName = unicode(magnitudeName)

    @classmethod
    # unit/factor mapping is defined for each class
    def magnitude(cls, name):
        """Returns a (numerical) magnitude matching a magnitude name"""
        name = unicode(name)
        try:
            return cls.magnitudeMapping[name]
        except KeyError:
            logging.warning(u"no matching magnitude to name {} found"
                    .format(name))

    @classproperty
    @classmethod
    # see _siMagnitudeName and thus cls.siMagnitudeName is defined for each subclass
    def siMagnitude(cls):
        return cls.magnitude(cls.siMagnitudeName)

    @classproperty
    @classmethod
    # see _siMagnitudeName is defined for each subclass
    def siMagnitudeName(cls):
        if cls._siMagnitudeName is None:
            raise NotImplementedError
        return cls._siMagnitudeName

    @property
    def displayMagnitude(self):
        return self.magnitude(self.displayMagnitudeName)

    @property
    # defined for instances only, given at init() time
    def displayMagnitudeName(self):
        return self._displayMagnitudeName

    @classproperty
    @classmethod
    # unit/factor mapping is defined for each class
    def magnitudeMapping(cls):
        if cls._magnitudeMap is None:
            raise NotImplementedError
        return cls._magnitudeMap

    @staticmethod
    # needs neither class not instance
    def invName(unitString):
        u""" Adds an ⁻¹ sign or removes it if already present"""
        unitString = unicode(unitString)
        if not u"⁻¹" in unitString:
            return unitString + u"⁻¹"
        else: 
            return unitString.replace( u"⁻¹", u"" )

    @property
    # for instances only, because of displayMagnitudeName
    def magnitudeConversion(self):
        """
        Scaling factor to move from display magnitude to si units.
        Required display argument:

            *displaymagnitudename* : The name of the magnitude to convert from

        Optional display argument:

            *simagnitudename* : The name of the magnitude to convert to.
                                Defaults to self.siMagnitudeName

        Returns:
            *float* : A scaling factor for display unit to scale to si unit.
        """
        # find display units:
        iUnit = self.magnitudeMapping[self.displayMagnitudeName]
        # find si units
        oUnit = self.magnitudeMapping[self.siMagnitudeName]
        return iUnit / oUnit

    def toSi(self, value):
        if isinstance(value, collections.Sequence): # for lists&tuples
            return type(value)((v * self.magnitudeConversion for v in value))
        # else:
        return value * self.magnitudeConversion

    def toDisplay(self, value):
        if isinstance(value, collections.Sequence): # for lists&tuples
            return type(value)((v / self.magnitudeConversion for v in value))
        # else:
        return value / self.magnitudeConversion

    @classmethod
    def name(cls):
        return cls.__name__

class Length(Unit):
    _magnitudeMap = {
        u"Å" : 1e-10,
        u"nm": 1e-9,
        u"µm": 1e-6,
        u"mm": 1e-3,
        u"cm": 1e-2,
        u"m" : 1e0
    }
    _siMagnitudeName = u"m"

NM = Length(u"nm")

class Area(Unit):
    _magnitudeMap = {
        u"Å²" : 1e-20,
        u"nm²": 1e-18,
        u"µm²": 1e-12,
        u"mm²": 1e-6,
        u"m²" : 1e0,
    }
    _siMagnitudeName = u"m²"

class Volume(Unit):
    _magnitudeMap = {
        u"Å³" : 1e-30,
        u"nm³": 1e-27,
        u"µm³": 1e-18,
        u"mm³": 1e-9,
        u"m³" : 1e0,
    }
    _siMagnitudeName = u"m³"

class Angle(Unit):
    _magnitudeMap = {
        u"°"  : pi/180.0, # unicode U+00B0
        u"'"  : pi/  3.0,
        u'"'  : pi/  0.05,
        u"rad": pi/  1.0,
    }
    _siMagnitudeName = u"rad"

class Temperature(Unit):
    _magnitudeMap = {
        u"K":  1.0,
    }
    _siMagnitudeName = u"K"

K = Temperature(u"K")

class Viscosity(Unit):
    _magnitudeMap = {
        u"mPa·s":  1e3,
        u"Pa·s":   1.0, # kg/(s·m)
    }
    _siMagnitudeName = u"Pa·s"

VIS = Viscosity(u"mPa·s")

class SLD(Unit):
    _magnitudeMap = {
        u"Å⁻²" : 1e20,
        u"nm⁻²": 1e18,
        u"µm⁻²": 1e12,
        u"mm⁻²": 1e6,
        u"cm⁻²": 1e4,
        u"m⁻²" : 1e0,
    }
    _siMagnitudeName = u"m⁻²"

class ScatteringVector(Unit):
    _magnitudeMap = {
        u"Å⁻¹" : 1e10,
        u"nm⁻¹": 1e9,
        u"µm⁻¹": 1e6,
        u"mm⁻¹": 1e3,
        u"cm⁻¹": 1e2,
        u"m⁻¹" : 1e0,
    }
    _siMagnitudeName = u"m⁻¹"

class ScatteringIntensity(Unit):
    _magnitudeMap = {
        u"(cm sr)⁻¹": 1e2,
        u"(m sr)⁻¹" : 1e0,
    }
    _siMagnitudeName = u"(m sr)⁻¹"

class Fraction(Unit):
    _magnitudeMap = {
        u"%": 1e-2,
        u"-": 1e0,
        u"" : 1e0,
    }
    _siMagnitudeName = u"-"

class NoUnit(Unit):
    _magnitudeMap = {
        u"" : 1e0,
        u"-": 1e0,
    }
    _siMagnitudeName = u"-"

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
# vim: set ts=4 sts=4 sw=4 tw=0:
