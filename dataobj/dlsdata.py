# -*- coding: utf-8 -*-
# dataobj/dlsdata.py

"""
Represents data from dynamic light scattering (DLS) measurement.

"""

from __future__ import absolute_import # PEP328

import logging
import numpy as np # For arrays
from utils import classproperty
from utils.units import Length, ScatteringVector, ScatteringIntensity, Angle

from dataobj.dataobj import DataObj

def _makeProperty(varName):
    def getter(selforcls):
        return getattr(selforcls, varName)
    return property(getter)

class DLSData(DataObj):
    """Represents one data set.
    """
    _properties = ("sampleName", "description", "correlation", "angles",
                   "temperature", "viscosity",
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

    @property
    def tau(self):
        return self.correlation[:, 0]

    def correlationAtAngle(self, angleIdx):
        return self.correlation[:, angleIdx]

    @property
    def count(self):
        return self.correlation.shape[0]

    # scattering angles

    def setAngles(self, angles):
        self._angles = angles

    @property
    def numAngles(self):
        return len(self.angles)

    def setTemperature(self, temp):
        self._temperature = temp

    def setViscosity(self, vis):
        self._viscosity = vis

    def setRefractiveIndex(self, refIdx):
        self._refractiveIndex = refIdx

    def setWaveLength(self, wavelen):
        self._wavelength = wavelen

    def __init__(self, **kwargs):
        super(DLSData, self).__init__(**kwargs)

    def __str__(self):
        out = [ "## " + self.__class__.__name__ + " ##" ]
        for p in self._properties:
            out.append(u"{0}: {1}".format(p, getattr(self, "_" + p)))
        return u"\n".join(out)

    @classmethod
    def setPropertyGetters(cls):
        for p in cls._properties:
            attr = "_" + p
            setattr(cls, attr, None) # init value
            setattr(cls, p, _makeProperty(attr))

DLSData.setPropertyGetters()

if __name__ == "__main__":
    import doctest
    doctest.testmod()

# vim: set ts=4 sts=4 sw=4 tw=0:
