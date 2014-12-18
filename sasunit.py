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

>>> rUnit = Length(simagnitudename = "m",
...                displaymagnitudename = "nm")
>>> rUnit.magnitudeConversion('nm', 'm')
1e-09

or:
>>> rUnit.magnitudeConversion('nm')
1e-09

Selecting a default: 
>>> qUnit = ScatteringVector(simagnitudename = u"m⁻¹",
...                          displaymagnitudename = u"cm⁻¹")
>>> qUnit.magnitudeConversion(u"cm⁻¹")
100.0

"""

import logging
import numpy as np # For arrays
from numpy import pi
from cutesnake.utils.tests import testfor
from cutesnake.utils.classproperty import classproperty

class SASUnit(object):
    _magnitudeDict = None
    _siMagnitudeName = u"" 
    _displayMagnitudeName = u""

    def __init__(self, **kwargs):
        """process display. Input should contain keywords defined above"""
        outArg = kwargs.get('simagnitudename', None)
        inArg = kwargs.get('displaymagnitudename', None)

        # TODO: assert given in/out args exist in selected magDict
#        import sys
#        print >>sys.__stderr__, "SASUnit:", inArg, outArg, self.invMagnitudeDict, outArg in self.invMagnitudeDict
#        print >>sys.__stderr__, self.invMagnitudeDict.items()
        if outArg is not None:
            testfor(outArg in self.invMagnitudeDict, ValueError,
                    "si magnitude name not in chosen magnitude dict!")
            self.siMagnitudeName = outArg
        if inArg is not None:
            testfor(inArg in self.invMagnitudeDict, ValueError,
                    "display magnitude name not in chosen magnitude dict!")
            self.displayMagnitudeName = inArg

    def magnitude(self, name):
        """Returns a (numerical) magnitude matching a magnitude name"""
        try:
            return self.invMagnitudeDict[name]
        except KeyError:
            logging.warning('no matching magnitude to name {} found'
                    .format(name))

    @property
    def siMagnitude(self):
        return self.magnitude(self.siMagnitudeName)

    @property
    def siMagnitudeName(self):
        return self._siMagnitudeName

    @siMagnitudeName.setter
    def siMagnitudeName(self, name):
        if name in self.invMagnitudeDict:
            self._siMagnitudeName = name
        else:
            logging.warning('no valid si magnitude name used.')

    @property
    def displayMagnitude(self):
        return self.magnitude(self.displayMagnitudeName)

    @property
    def displayMagnitudeName(self):
        return self._displayMagnitudeName

    @displayMagnitudeName.setter
    def displayMagnitudeName(self, name):
        if name in self.invMagnitudeDict:
            self._displayMagnitudeName = name
        else:
            logging.warning('no valid display magnitude name: {}'.format(name))

    @classproperty
    @classmethod
    def invMagnitudeDict(cls):
        return cls._magnitudeDict

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
        iUnit = self.invMagnitudeDict[displaymagnitudename.strip()]
        #find si units
        oUnit = self.invMagnitudeDict[simagnitudename.strip()]
        return iUnit / oUnit

class Length(SASUnit):
    _magnitudeDict = {
        u"Å" : 1e-10,
        u"nm": 1e-9,
        u"µm": 1e-6,
        u"mm": 1e-3,
        u"cm": 1e-2,
        u"m" : 1e0
    }

class Area(SASUnit):
    _magnitudeDict = {
        u"Å²" : 1e-20,
        u"nm²": 1e-18,
        u"µm²": 1e-12,
        u"mm²": 1e-6,
        u"m²" : 1e0,
    }

class Volume(SASUnit):
    _magnitudeDict = {
        u"Å³" : 1e-30,
        u"nm³": 1e-27,
        u"µm³": 1e-18,
        u"mm³": 1e-9,
        u"m³" : 1e0,
    }

class Angle(SASUnit):
    _magnitudeDict = {
        u"˚"  : 180.0/pi,
        u"'"  :   3.0/pi,
        u'"'  :   0.05/pi,
        u"rad":   1.0,
    }

class SLD(SASUnit):
    _magnitudeDict = {
        u"Å⁻²" : 1e20,
        u"nm⁻²": 1e18,
        u"µm⁻²": 1e12,
        u"mm⁻²": 1e6,
        u"cm⁻²": 1e4,
        u"m⁻²" : 1e0,
    }

class ScatteringVector(SASUnit):
    _magnitudeDict = {
        u"Å⁻¹" : 1e10,
        u"nm⁻¹": 1e9,
        u"µm⁻¹": 1e6,
        u"mm⁻¹": 1e3,
        u"cm⁻¹": 1e2,
        u"m⁻¹" : 1e0,
    }

class ScatteringIntensity(SASUnit):
    _magnitudeDict = {
        u"(cm sr)⁻¹": 1e2,
        u"(m sr)⁻¹" : 1e0,
    }

class Fraction(SASUnit):
    _magnitudeDict = {
        u"%": 1e-2,
        u"-": 1e0,
        u"" : 1e0,
    }

class NoUnit(SASUnit):
    _magnitudeDict = {
        u"" : 1e0,
        u"-": 1e0,
    }

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
# vim: set ts=4 sts=4 sw=4 tw=0:
