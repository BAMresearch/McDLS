# -*- coding: utf-8 -*-
# dataobj/dlsdata.py

"""
Represents data from dynamic light scattering (DLS) measurement.

"""

from __future__ import absolute_import # PEP328

import logging
import copy
from collections import OrderedDict
from numpy import (pi, sin, array, dstack, hstack, newaxis, repeat, outer,
                   flipud, concatenate, empty)
from utils import classproperty, isCallable, isInteger, isList
from utils.units import (Length, ScatteringVector, ScatteringIntensity, Angle)
from dataobj import DataObj, DataVector, DataConfig
from models import DLSModel

# Boltzmann constant in m²·kg·s⁻²·K⁻¹ (SI units)
KB = 1.38064852 * 1e-23

def _makeProperty(varName):
    def getter(selforcls):
        return getattr(selforcls, varName)
    return property(getter)

def _privPropName(propName):
    return "_" + propName

def _propSetterName(propName):
    return "set" + propName[0].upper() + propName[1:]

class DLSConfig(DataConfig):
    def __init__(self):
        super(DLSConfig, self).__init__()
        self.nBin.setValue(0)

    def overrideDefaults(self):
        # clip preceding 2 points by default
        self.x0LowClip.setValue(2)

    @property
    def showParams(self):
        lst = super(DLSConfig, self).showParams
        lst.remove("nBin")
        return lst

DLSConfig.factory()

class MultiDataVector(DataVector):
    _count = None
    _wasRepeated = False

    def __init__(self, name, raw, rawU = None, count = None, **kwargs):
        # not beautiful but the easiest for now ...
        super(MultiDataVector, self).__init__(name,
                self.flatten(raw, count),
                self.flatten(rawU, count), **kwargs)

    # http://stackoverflow.com/a/7020271
    def flatten(self, a, count = None):
        if a is None:
            return a
        if (a.ndim == 1 or min(a.shape) == 1):
            a = a.reshape((-1, 1))
            if isInteger(count):
                # repeat one-dimensional data by <count> columns
                a = repeat(a, count, axis = 1)
                self._wasRepeated = True
        self._count = a.shape[1]
        a = self.concat(a) # concat the columns behind each other
        return a

    @property
    def rawDataSrcShape(self):
        return self.unflatten(self.rawData)

    @property
    def rawDataUSrcShape(self):
        return self.unflatten(self.rawDataU)

    @property
    def siDataSrcShape(self):
        return self.unflatten(self.siData)

    def unflatten(self, a):
        if a is None:
            return a
        if self._count == 1:
            return a.reshape((-1, 1)) # nothing to do
        rowCount = len(a) / self._count
        if self._wasRepeated: # just take the first non-duplicates
            result = a[:rowCount]
        else:
            result = self.unconcat(a, rowCount, self._count)
        return result

    @classmethod
    def concat(cls, a):
        """Returns a flat array, column-wise concatenated and each second
        column reversed. This way the data stays continuous and the plotting
        goes back and forth through the domain of definition as often as there
        are angles."""
        if a.shape[1] == 1:
            return a.flatten() # nothing to do
        o = a.copy()
        o[:,1::2] = o[::-1,1::2] # reverses the odd columns, see numpy.flipud()
        return o.ravel(order = 'F') # column-wise concat, Fortran style

    @classmethod
    def unconcat(cls, a, rows, cols):
        o = a.reshape((rows, cols), order = 'F').copy()
        o[:,1::2] = o[::-1,1::2] # reverses the odd columns, see numpy.flipud()
        return o

class DLSData(DataObj):
    """Represents one data set.
    """
    _properties = ("sampleName", "description",
                   "angles", "anglesUnit", "temperature", "viscosity",
                   "refractiveIndex", "wavelength",
                   "measIndices",
                   # calculated properties
                   "scatteringVector", "gammaDivR", "tauGamma")

    @classproperty
    @classmethod
    def displayDataDescr(cls):
        return ("Title", "Measurements", "Data points",
                "# Angles", "Angle(s)", "Description")

    @classproperty
    @classmethod
    def displayData(cls):
        return ("title", "measIndicesStr", "count",
                "numAngles", "anglesToStr", "description")

    @property
    def dataContent(self):
        """Shows the content of the loaded data: DLS?"""
        return ""

    @classproperty
    @classmethod
    def sourceName(cls):
        return "Dynamic Light Scattering"

    @property
    def configType(self):
        return DLSConfig

    @property
    def modelType(self):
        return DLSModel

    # object specific properties

    def setSampleName(self, sampleName):
        self._sampleName = sampleName

    def setSampleDescription(self, descr):
        self._description = descr

    def setMeasIndices(self, measIndices):
        """Sets the measurement index of this data set. Expects a list of
        tuples containing two integers of the measurement group index and the
        measurement index within that group."""
        def verifyMeasIndex(measIndex):
            assert isList(measIndex) and len(measIndex) == 2
            assert all([isInteger(i) for i in measIndex])
        assert isList(measIndices)
        [verifyMeasIndex(mi) for mi in measIndices]
        self._measIndices = measIndices

    @property
    def measIndicesStr(self):
        summary = OrderedDict()
        for g, i in self.measIndices:
            if g not in summary:
                summary[g] = []
            summary[g].append(i)
        res = ";".join(["{group}: {indices}".format(group = g,
            indices = ",".join([str(i) for i in lst]))
                for g, lst in summary.iteritems()])
        return res

    # define DataObj interface

    @property
    def tau(self):
        return self.x0

    @property
    def correlation(self):
        return self.f

    def setTau(self, tauUnit, rawArray):
        self.x0 = MultiDataVector(u"τ", rawArray.flatten(), unit = tauUnit,
                                  count = self.numAngles)
        self._calcTauGamma()

    def setCorrelation(self, rawArray, rawArrayU = None):
        assert self.isValidInput(rawArray), "Invalid data from file!"
        assert rawArray.shape[1] == self.numAngles, (
            "Correlation intensity: {} columns != {} scattering angles"
            .format(rawArray.shape[1], self.numAngles))
        if rawArrayU is not None:
            assert self.isValidInput(rawArrayU), "Invalid data from file!"
            assert rawArray.shape == rawArrayU.shape, \
                "Shapes of the provided data and its uncertainties do not match!"
        self.f = MultiDataVector(u"G_2(τ)-1", rawArray, rawU = rawArrayU)

    @property
    def count(self):
        return len(self.x0.sanitized)

    # scattering angles

    def setAngles(self, anglesUnit, angles):
        """Expects angles in si units and their associated unit for proper
        formatting."""
        self._angles = angles
        self._anglesUnit = anglesUnit
        self._calcScatteringVector()

    @property
    def numAngles(self):
        if self.angles is None:
            return 0
        return len(self.angles)

    @property
    def anglesToStr(self):
        return self.anglesFmt()

    def anglesFmt(self, fmt = None):
        unit = self.anglesUnit
        if fmt is None:
            fmt = u"{0:.0f}{1}"
        return u";".join((
            unicode(fmt).format(unit.toDisplay(a), unit.displayMagnitudeName)
                    for a in self.angles))

    # temperature, viscosity, refractiveIndex, wavelength including std.err

    def setTemperature(self, temp, stddev = None):
        self._temperature = (temp, stddev)
        self._calcGammaDivR()

    def setViscosity(self, vis, stddev = None):
        self._viscosity = (vis, stddev)
        self._calcGammaDivR()

    def setRefractiveIndex(self, refIdx, stddev = None):
        self._refractiveIndex = (refIdx, stddev)
        self._calcScatteringVector()

    def setWavelength(self, wavelen, stddev = None):
        self._wavelength = (wavelen, stddev)
        self._calcScatteringVector()

    # calculated properties, prepared for model evaluation

    def _calcScatteringVector(self):
        """Calculates the scattering vector which is part of the scattering
        formula."""
        if (self.refractiveIndex is None or
            self.angles is None or
            self.wavelength is None):
            return
        self._scatteringVector = (
            4. * pi * self.refractiveIndex[0] * sin(self.angles * .5)
                / self.wavelength[0])
        self._calcGammaDivR()

    def _calcGammaDivR(self):
        """Calculates the gamma value without considering the hydrodynamical
        radius which is contributed (divided) by the model later on."""
        if (self.temperature is None or
            self.viscosity is None or
            self.scatteringVector is None):
            return
        q2 = self._scatteringVector * self._scatteringVector
        # 0.5 ensures the necessary sqrt() within the model
        self._gammaDivR = ( (- q2 * self.temperature[0] * KB)
                         / (6. * pi * self.viscosity[0]) )
        self._calcTauGamma()

    def _calcTauGamma(self):
        """tau premultiplied by gamma without hydrodynamical radius produces
        a matrix of the same dimensions as the correlation data."""
        if self.tau is None or self.gammaDivR is None:
            return
        self._tauGammaMat = outer(self.tau.siDataSrcShape, self.gammaDivR)
        self._tauGamma = MultiDataVector(u"tauGamma", self._tauGammaMat)

    def _propagateMask(self, *args):
        super(DLSData, self)._propagateMask(*args)
        self._tauGamma.validIndices = self.tau.validIndices

    def accumulate(self, others):
        """Combines several measured data of the same sample into a single
        data set with an uncertainty for each measured value."""
        # consider only data of the same type and sample name
        others = [o for o in others if 
                    o is not None and
                    isinstance(o, DLSData) and
                    o.sampleName == self.sampleName]
        if not len(others):
            return None
        if len(others) == 1:
            return self
        # average basic properties
        for prop in ("temperature", "viscosity",
                     "refractiveIndex", "wavelength"):
            arr = array([getattr(o, prop)[0] for o in others])
            setFunc = getattr(self, _propSetterName(prop))
            if isCallable(setFunc):
                setFunc(arr.mean(), arr.std())
        # combine all measurement indices
        self.setMeasIndices(tuple((mi for o in others for mi in o.measIndices)))
        # angles unchanged, but ensure they're identical everywhere
        assert array([self.angles == o.angles for o in others]).all(), \
               "Scattering angles differ between all DLS data to be combined!"
        # calculate average correlation values and their standard deviation
        stacked = dstack((o.correlation.rawDataSrcShape for o in others))
        # combine the mean across all data sets with the existing tau
        # combine the std. deviation with the existing tau
        self.setCorrelation(stacked.mean(-1), stacked.std(-1))
        assert len(self.x0.siData) == len(self.f.siData), \
            "Dimensions of flattened data arrays do not match!"
        # reset config in order to fix callbacks
        self.setConfig(self.configType())
        return self

    def splitPerAngle(self):
        lst = []
        if self.numAngles == 1:
            return [self]
        for i in range(self.numAngles):
            another = copy.copy(self)
            another.setAngles(self.anglesUnit, self.angles[i, newaxis])
            another.setTau(self.tau.unit, self.tau.rawDataSrcShape)
            another.setCorrelation(
                    self.correlation.rawDataSrcShape[:, i, newaxis],
                    self.correlation.rawDataUSrcShape[:, i, newaxis])
            # reset config in order to fix callbacks
            another.setConfig(another.configType())
            lst.append(another)
        return lst

    def ids(self):
        """Returns a text containing the objects ids of the embedded data
        objects. (for debugging purposes)"""
        out = [repr(self)]
        for p in self._properties:
            out.append("{}: {}".format(p, id(getattr(self, _privPropName(p)))))
        return u"\n".join(out)

    def __init__(self, **kwargs):
        super(DLSData, self).__init__(**kwargs)

    def setConfig(self, config):
        if not super(DLSData, self).setConfig(config):
            return # no update, nothing todo
        self.config.x0Low.setUnit(self.tau.unit)
        self.config.x0High.setUnit(self.tau.unit)

    def __str__(self):
        out = [u"## {0} '{1}'".format(self.__class__.__name__, self.title)]
        for p in self._properties:
            out.append(u"{0}: {1}".format(p, getattr(self, _privPropName(p))))
        return u"\n".join(out)

    def __hash__(self):
        value = hash(self.title) ^ hash(self.filename)
        for p in self._properties:
            propData = getattr(self, _privPropName(p))
            try:
                value ^= hash(propData)
            except TypeError: # numpy.ndarray
                value ^= hash(propData.data.tobytes())
        return value

    @classmethod
    def setPropertyGetters(cls):
        for p in cls._properties:
            attr = _privPropName(p)
            setattr(cls, attr, None) # init value = None
            setattr(cls, p, _makeProperty(attr))

DLSData.setPropertyGetters()

if __name__ == "__main__":
    import doctest
    doctest.testmod()

# vim: set ts=4 sts=4 sw=4 tw=0:
