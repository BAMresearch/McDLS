# -*- coding: utf-8 -*-
# mcsas/sasdata.py

"""
Represents data associated with a measurement by small angle scattering (SAS).
Some examples and tests.

>>> import numpy
>>> testdata = numpy.random.rand(4,4)
>>> testtitle = "some title"
>>> from sasdata import SASData

Testing copy()
>>> first = SASData(testtitle, testdata)
>>> first.title == testtitle
True
>>> numpy.all(first.origin == testdata)
True
>>> second = first.copy()
>>> second.title == testtitle
True
>>> numpy.all(second.origin == testdata)
True
>>> first == second
True
"""

import logging
import numpy as np # For arrays

class SASUnits(object):
    """
    Defines methods for using and manipulating units of variables. 
    Some default magnitude-name dictionaries are provided, but the user can 
    supply their own dictionary if required. Default unit to translate to 
    must be set.
    Required keyword arguments:
    *magnitudedict*: a dictionary of magnitude - name pairs. Names must be 
        unicode strings.
    *outputmagnitudename*: the output magnitude name

    Available methods:

    """
    _magnitudeDict = dict()
    _outputMagnitudeName = u"" 
    _inputMagnitudeName = u""
    #default library
    _lengthMagnitudeDict = { 1e-10 : u"Å",
            1e-9 : u"nm",
            1e-6 : u"µm",
            1e-3 : u"nm",
            1e-2 : u"cm",
            1e0 : u"m"}

    def __init__(self, **kwargs):
        """process input. Input should contain keywords defined above"""
        dictArg = kwargs.get('magnitudedict', None)
        outArg = kwargs.get('outputmagnitudename', None)
        inArg = kwargs.get('inputmagnitudename', None)

        if dictArg is not None:
            self.magnitudeDict = dictArg
        if outArg is not None:
            self.outputMagnitudeName = outArg
        if inArg is not None:
            self.inputMagnitudeName = inArg

    def magnitude(self, name):
        """Returns a (numerical) magnitude matching a magnitude name"""
        try:
            return self.invMagnitudeDict()[name]
        except KeyError:
            logging.warning('no matching magnitude to name {} found'
                    .format(name))

    @property
    def outputMagnitude(self):
        return self.magnitude(self.outputMagnitudeName)
    @outputMagnitude.setter
    def outputMagnitude(self, magnum):
        try:
            self.outputMagnitudeName = self.magnitudeDict[magnum]
        except KeyError:
            logging.warning('no matching magnitude name for output magnitude {}'
                    .format(magname))

    @property
    def outputMagnitudeName(self):
        return self._outputMagnitudeName
    @outputMagnitudeName.setter
    def outputMagnitudeName(self, name):
        if name in self.invMagnitudeDict():
            self._outputMagnitudeName = name
        else:
            logging.warning('no valid output magnitude name used.')

    @property
    def inputMagnitude(self):
        return self.magnitude(self.inputMagnitudeName)
    @inputMagnitude.setter
    def inputMagnitude(self, magnum):
        try:
            self.inputMagnitudeName = self.magnitudeDict[magnum]
        except KeyError:
            logging.warning('no matching magnitude name for input magnitude {}'
                    .format(magname))

    @property
    def inputMagnitudeName(self):
        return self._inputMagnitudeName
    @inputMagnitudeName.setter
    def inputMagnitudeName(self, name):
        if name in self.invMagnitudeDict():
            self._inputMagnitudeName = name
        else:
            logging.warning('no valid input magnitude name used.')

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

    def magnitudeConversion(self, inputmagnitudename = None, 
            outputmagnitudename = None):
        """
        Scaling factor to move from input magnitude to output units.
        Required input argument:
        *inputmagnitudename* : The name of the magnitude to convert from
        Optional input argument:
        *outputmagnitudename* : The name of the magnitude to convert to. 
            Defaults to self.outputMagnitudeName

        Returns: 
        *float* : A scaling factor for input unit to scale to output unit.
        """
        if outputmagnitudename is None:
            outputmagnitudename = self.outputMagnitudeName

        #find input units:
        iUnit = self.invMagnitudeDict()[inputmagnitudename.strip()]
        #find output units
        oUnit = self.invMagnitudeDict()[outputmagnitudename.strip()]
        return iUnit / oUnit
    
# vim: set ts=4 sts=4 sw=4 tw=0:
