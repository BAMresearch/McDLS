# -*- coding: utf-8 -*-
# dataobj/dlsdata.py

"""
Represents data from dynamic light scattering (DLS) measurement.

"""

from __future__ import absolute_import # PEP328

import logging
from numpy import array as np_array
from numpy import dstack as np_dstack
from numpy import hstack as np_hstack
from numpy import newaxis as np_newaxis
from utils import classproperty, isFunction
from utils.units import Length, ScatteringVector, ScatteringIntensity, Angle
from dataobj.dataobj import DataObj

import sys

def _makeProperty(varName):
    def getter(selforcls):
        return getattr(selforcls, varName)
    return property(getter)

def _privPropName(propName):
    return "_" + propName

def _propSetterName(propName):
    return "set" + propName[0].upper() + propName[1:]

class DLSData(DataObj):
    """Represents one data set.
    """
    _properties = ("sampleName", "description",
                   "correlation", "correlationError",
                   "angles", "temperature", "viscosity",
                   "refractiveIndex", "wavelength",)

    @classproperty
    @classmethod
    def displayDataDescr(cls):
        return ("Filename ", "Data points ", "Angles ",
                "Sample Name", "Description")

    @classproperty
    @classmethod
    def displayData(cls):
        return ("title", "count", "numAngles", "sampleName", "description")

    @property
    def dataContent(self):
        """Shows the content of the loaded data: DLS?"""
        return ""

    def setSampleName(self, sampleName):
        self._sampleName = sampleName

    def setSampleDescription(self, descr):
        self._description = descr

    # correlation data

    def setCorrelation(self, rawArray):
        assert self.isValidInput(rawArray), "Invalid data from file!"
        self._correlation = rawArray

    def setCorrelationError(self, rawArray):
        assert self.isValidInput(rawArray), "Invalid data from file!"
        self._correlationError = rawArray

    @property
    def tau(self):
        return self.correlation[:, 0]

    def correlationAtAngle(self, angleIdx):
        return self.correlation[:, angleIdx + 1]

    @property
    def count(self):
        return self.correlation.shape[0]

    # scattering angles

    def setAngles(self, angles):
        self._angles = angles

    @property
    def numAngles(self):
        return len(self.angles)

    # temperature, viscosity, refractiveIndex, wavelength including std.err

    def setTemperature(self, temp, stddev = None):
        self._temperature = (temp, stddev)

    def setViscosity(self, vis, stddev = None):
        self._viscosity = (vis, stddev)

    def setRefractiveIndex(self, refIdx, stddev = None):
        self._refractiveIndex = (refIdx, stddev)

    def setWavelength(self, wavelen, stddev = None):
        self._wavelength = (wavelen, stddev)

    def accumulate(self, others):
        others = [o for o in others if 
                    o is not None and
                    isinstance(o, DLSData) and
                    o.sampleName == self.sampleName]
        self.title = self.sampleName + " (averaged)"
        for prop in ("temperature", "viscosity",
                     "refractiveIndex", "wavelength"):
            array = np_array([getattr(o, prop)[0] for o in others])
            setFunc = getattr(self, _propSetterName(prop))
            if isFunction(setFunc):
                setFunc(array.mean(), array.std())
        assert np_array([self.angles == o.angles for o in others]).all(), \
               "Scattering angles differ between all DLS data to be combined!"
        # calculate average correlation values and their standard deviation
        stacked = np_dstack((o.correlation for o in others))
        # combine the mean across all data sets with the existing tau
        corr = np_hstack((self.tau[:, np_newaxis], stacked.mean(2)[:, 1:]))
        # combine the std. deviation with the existing tau
        stddev = np_hstack((self.tau[:, np_newaxis], stacked.std(2)[:, 1:]))
        self.setCorrelation(corr)
        self.setCorrelationError(stddev)
#        print >>sys.__stderr__, unicode(self)
        return self

    def __init__(self, **kwargs):
        super(DLSData, self).__init__(**kwargs)

    def __str__(self):
        out = []
        out.append(u"## {0} '{1}'".format(self.__class__.__name__, self.title))
        for p in self._properties:
            out.append(u"{0}: {1}".format(p, getattr(self, _privPropName(p))))
        return u"\n".join(out)

    @classmethod
    def setPropertyGetters(cls):
        for p in cls._properties:
            attr = _privPropName(p)
            setattr(cls, attr, None) # init value
            setattr(cls, p, _makeProperty(attr))

DLSData.setPropertyGetters()

if __name__ == "__main__":
    import doctest
    doctest.testmod()

# vim: set ts=4 sts=4 sw=4 tw=0:
