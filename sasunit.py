# -*- coding: utf-8 -*-
# mcsas/sasdata.py

"""
Defines methods for using and manipulating units of variables. 
Some default magnitude-name dictionaries are provided, but the user can 
supply their own dictionary if required. Default unit to translate to 
must be set.
Required keyword arguments:
*magnitudedict*: a dictionary of magnitude - name pairs. Names must be 
    unicode strings.
*simagnitudename*: the si magnitude name

Example usage: 

>>> rUnit = SASUnit(magnitudedict = {1e-9 : u"nm", 1e0 : u"m"}, 
    simagnitudename = "m", 
    displaymagnitudename = "nm")
>>> rUnit.magnitudeConversion('nm', 'm')
or:
>>> rUnit.magnitudeConversion('nm')
1e-9

Selecting a default: 
>>> qUnit = SASUnit(magnitudedict = 'q', simagnitudename = u"m⁻¹", displaymagnitudename = u"cm⁻¹")
>>> qUnit.magnitudeConversion(u"cm⁻¹")
100.0

"""

import logging
import numpy as np # For arrays

class SASUnit(object):
    _magnitudeDict = dict()
    _siMagnitudeName = u"" 
    _displayMagnitudeName = u""
    #default library, if growing out of bounds should be put in json dict
    _defaultDicts = {
            'length' : { 1e-10 : u"Å", 
                1e-9 : u"nm",
                1e-6 : u"µm",
                1e-3 : u"mm",
                1e-2 : u"cm",
                1e0 : u"m"},
            'area' : { 1e-20 : u"Å²",
                1e-18 : u"nm²",
                1e-12 : u"µm²",
                 1e-6 : u"mm²",
                 1.   : u"m²"}
            'volume' : { 1e-30 : u"Å³",
                1e-27 : u"nm³",
                1e-18 : u"µm³",
                 1e-9 : u"mm³",
                   1. : u"m³"}
            'SLD' : { 1e20 : u"Å⁻²",
                1e18 : u"nm⁻²",
                1e12 : u"µm⁻²",
                 1e6 : u"mm⁻²",
                 1e4 : u"cm⁻²",
                  1. : u"m⁻²"}
            'q' :  { 1e10 : u"Å⁻¹",
                1e9 : u"nm⁻¹",
                1e6 : u"µm⁻¹",
                1e3 : u"mm⁻¹",
                1e2 : u"cm⁻¹",
                1e0 : u"m⁻¹"},
            'I' :  { 1e2 : u"(cm sr)⁻¹",
                1e0 : u"(m sr)⁻¹"},
            'none' : { 1. : u"",
                1. : u"-"}

            }

    def __init__(self, **kwargs):
        """process display. Input should contain keywords defined above"""
        dictArg = kwargs.get('magnitudedict', None)
        outArg = kwargs.get('simagnitudename', None)
        inArg = kwargs.get('displaymagnitudename', None)

        if dictArg is not None:
            if isinstance(dictArg, dict):
                #set the dict
                self.magnitudeDict = dictArg
            else: 
                #select one of the predefined dicts
                self.magnitudeDict = self.defaultDict(dictArg)
        if outArg is not None:
            self.siMagnitudeName = outArg
        if inArg is not None:
            self.displayMagnitudeName = inArg

    def defaultDict(self, dictArg = None):
        if dictArg in self._defaultDicts.keys():
            return self._defaultDicts[dictArg]
        else:
            logging.warning('dictionary for {} not found'.format(dictArg))
        logging.info('Default dictionaries available for: {}'
                .format(self._defaultDicts.keys()))
        return None

    def magnitude(self, name):
        """Returns a (numerical) magnitude matching a magnitude name"""
        try:
            return self.invMagnitudeDict()[name]
        except KeyError:
            logging.warning('no matching magnitude to name {} found'
                    .format(name))

    @property
    def siMagnitude(self):
        return self.magnitude(self.siMagnitudeName)
    @siMagnitude.setter
    def siMagnitude(self, magnum):
        try:
            self.siMagnitudeName = self.magnitudeDict[magnum]
        except KeyError:
            logging.warning('no matching magnitude name for si magnitude {}'
                    .format(magname))

    @property
    def siMagnitudeName(self):
        return self._siMagnitudeName
    @siMagnitudeName.setter
    def siMagnitudeName(self, name):
        if name in self.invMagnitudeDict():
            self._siMagnitudeName = name
        else:
            logging.warning('no valid si magnitude name used.')

    @property
    def displayMagnitude(self):
        return self.magnitude(self.displayMagnitudeName)
    @displayMagnitude.setter
    def displayMagnitude(self, magnum):
        try:
            self.displayMagnitudeName = self.magnitudeDict[magnum]
        except KeyError:
            logging.warning('no matching magnitude name for display magnitude {}'
                    .format(magname))

    @property
    def displayMagnitudeName(self):
        return self._displayMagnitudeName
    @displayMagnitudeName.setter
    def displayMagnitudeName(self, name):
        if name in self.invMagnitudeDict():
            self._displayMagnitudeName = name
        else:
            logging.warning('no valid display magnitude name used.')

    @property
    def magnitudeDict(self):
        return self._magnitudeDict
    @magnitudeDict.setter
    def magnitudeDict(self, dictionary):
        if isinstance(dictionary, dict):
            self._magnitudeDict = dictionary
        #defaults not working yet
        #elif isStr(dictionary):
        #    if dictionary.lower() == 'length':
        #        self._MagnitudeDict = self._lengthMagnitudeDict

    def invMagnitudeDict(self):
        return {v:k for k, v in self.magnitudeDict.items()}

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
        iUnit = self.invMagnitudeDict()[displaymagnitudename.strip()]
        #find si units
        oUnit = self.invMagnitudeDict()[simagnitudename.strip()]
        return iUnit / oUnit
    
# vim: set ts=4 sts=4 sw=4 tw=0:
