# -*- coding: utf-8 -*-
# mcsas/sasdata.py

u"""
Defines methods for using and manipulating units of variables. 
Some default magnitude-name dictionaries are provided, but the user can 
supply their own dictionary if required. Default unit to translate to 
must be set.
Required keyword arguments:
*magnitudedict*: a dictionary of magnitude - name pairs. Names must be 
    unicode strings.
*simagnitudename*: the si magnitude name

Example usage: 

>>> rUnit = Length("nm")
>>> rUnit.siMagnitudeName
u'm'
>>> rUnit.displayMagnitudeName
u'nm'
>>> rUnit.magnitudeConversion('nm', 'm')
1e-09

or:
>>> rUnit.magnitudeConversion('nm')
1e-09

Selecting a default: 
>>> qUnit = ScatteringVector(u"cm⁻¹")
>>> qUnit.magnitudeConversion(u"cm⁻¹")
100.0

"""

import logging
import numpy as np # For arrays
from numpy import pi
from cutesnake.utils.tests import testfor
from cutesnake.utils.classproperty import classproperty

class SASUnit(object):
    _magnitudeMap = None
    _siMagnitudeName = None
    _displayMagnitudeName = None

    def __init__(self, magnitudeName = None):
        """Set up a unit of measurement. The provided magnitude name is used
        by default in the UI. Using SI units if unspecified."""
        if magnitudeName is None:
            magnitudeName = self.siMagnitudeName
        elif magnitudeName not in self.magnitudeMapping:
            logging.warning("Provided default magnitude name '{mn}' is not "
                            "available: '{map}'!".format(
                            nm = magnitudeName, map = self.magnitudeMapping))
        self._displayMagnitudeName = unicode(magnitudeName)

    def magnitude(self, name):
        """Returns a (numerical) magnitude matching a magnitude name"""
        try:
            return self.magnitudeMapping[name]
        except KeyError:
            logging.warning('no matching magnitude to name {} found'
                    .format(name))

    @property
    def siMagnitude(self):
        return self.magnitude(self.siMagnitudeName)

    @property
    def siMagnitudeName(self):
        if self._siMagnitudeName is None:
            raise NotImplementedError
        return self._siMagnitudeName

    @property
    def displayMagnitude(self):
        return self.magnitude(self.displayMagnitudeName)

    @property
    def displayMagnitudeName(self):
        return self._displayMagnitudeName

    @classproperty
    @classmethod
    def magnitudeMapping(cls):
        if cls._magnitudeMap is None:
            raise NotImplementedError
        return cls._magnitudeMap

    def invName(self, unitString):
        u""" Adds an ⁻¹ sign or removes it if already present"""
        if not u"⁻¹" in unitString:
            return unitString + u"⁻¹"
        else: 
            return unitString.replace( u"⁻¹", u"" )

    def magnitudeConversion(self, displaymagnitudename = None, 
            simagnitudename = None):
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
        if displaymagnitudename is None:
            displaymagnitudename = self.displayMagnitudeName
        if simagnitudename is None:
            simagnitudename = self.siMagnitudeName

        #find display units:
        iUnit = self.magnitudeMapping[displaymagnitudename.strip()]
        #find si units
        oUnit = self.magnitudeMapping[simagnitudename.strip()]
        return iUnit / oUnit

class Length(SASUnit):
    _magnitudeMap = {
        u"Å" : 1e-10,
        u"nm": 1e-9,
        u"µm": 1e-6,
        u"mm": 1e-3,
        u"cm": 1e-2,
        u"m" : 1e0
    }
    _siMagnitudeName = u"m"

class Area(SASUnit):
    _magnitudeMap = {
        u"Å²" : 1e-20,
        u"nm²": 1e-18,
        u"µm²": 1e-12,
        u"mm²": 1e-6,
        u"m²" : 1e0,
    }
    _siMagnitudeName = u"m²"

class Volume(SASUnit):
    _magnitudeMap = {
        u"Å³" : 1e-30,
        u"nm³": 1e-27,
        u"µm³": 1e-18,
        u"mm³": 1e-9,
        u"m³" : 1e0,
    }
    _siMagnitudeName = u"m³"

class Angle(SASUnit):
    _magnitudeMap = {
        u"˚"  : 180.0/pi,
        u"'"  :   3.0/pi,
        u'"'  :   0.05/pi,
        u"rad":   1.0,
    }
    _siMagnitudeName = u"rad"

class SLD(SASUnit):
    _magnitudeMap = {
        u"Å⁻²" : 1e20,
        u"nm⁻²": 1e18,
        u"µm⁻²": 1e12,
        u"mm⁻²": 1e6,
        u"cm⁻²": 1e4,
        u"m⁻²" : 1e0,
    }
    _siMagnitudeName = u"m⁻²"

class ScatteringVector(SASUnit):
    _magnitudeMap = {
        u"Å⁻¹" : 1e10,
        u"nm⁻¹": 1e9,
        u"µm⁻¹": 1e6,
        u"mm⁻¹": 1e3,
        u"cm⁻¹": 1e2,
        u"m⁻¹" : 1e0,
    }
    _siMagnitudeName = u"m⁻¹"

class ScatteringIntensity(SASUnit):
    _magnitudeMap = {
        u"(cm sr)⁻¹": 1e2,
        u"(m sr)⁻¹" : 1e0,
    }
    _siMagnitudeName = u"(m sr)⁻¹"

class Fraction(SASUnit):
    _magnitudeMap = {
        u"%": 1e-2,
        u"-": 1e0,
        u"" : 1e0,
    }
    _siMagnitudeName = u"-"

class NoUnit(SASUnit):
    _magnitudeMap = {
        u"" : 1e0,
        u"-": 1e0,
    }
    _siMagnitudeName = u"-"

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
# vim: set ts=4 sts=4 sw=4 tw=0:
