# -*- coding: utf-8 -*-
# dataobj/dlsdata.py

"""
Represents data from dynamic light scattering (DLS) measurement.

"""

import logging
import copy
from collections import OrderedDict
from itertools import groupby, chain
from operator import itemgetter
import numpy as np
from numpy import (pi, sin, array, dstack, hstack, newaxis, repeat, outer,
                   flipud, concatenate, empty, zeros_like)

from ..utils import classproperty, isCallable, isInteger, isList, isIterable, hashNumpyArray
from ..utils.units import (Length, ScatteringVector, ScatteringIntensity, Angle,
                         NoUnit)
from ..bases.algorithm import Parameter
from . import DataObj, DataVector, DataConfig
from ..bases.model import DLSModel

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
    shortName = "DLS data configuration"
    parameters = (
        Parameter("plotCountRate", True, unit = NoUnit(),
            displayName = "plot count rate",
            description = "Plots the mean count rate with an error bar for the "
                "standard deviation of multiple measurements. "
                "It appears behind the correlation plot."),
        Parameter("doAverage", True, unit = NoUnit(),
            displayName = "average each angle",
            description = "When loading DLS data files, average the "
                "correlation data of each angle if there are multiple "
                "measurements of the same sample. <br />Finally, "
                "there will be a single data object for each angle which "
                "contains the correlation mean and its standard deviation "
                "interpreted as measurement uncertainty."),
        Parameter("outlierThreshold", 3.5, unit = NoUnit(),
            displayName = "outlier threshold",
            valueRange = (0., 1e3),
            description = "Modified Z-score above which measurements are "
                "rejected as outliers by using the median absolute deviation "
                "of all measurements for the same angle."),
    )
    medianCountRate = None
    # filterMask and outlierMap are complementary
    filterMask = None # for selecting good data during accumulate()
    outlierMap = None # for marking bad data in plotting

    def __init__(self):
        super(DLSConfig, self).__init__()
        self.nBin.setValue(0)

    @property
    def showParams(self):
        lst = super(DLSConfig, self).showParams
        lst.remove("nBin")
        return lst

    def update(self, other):
        """Override AlgorithmBase.update() to copy DataConfig values as well.
        """
        super(DLSConfig, self).update(other)
        self.medianCountRate = other.medianCountRate
        self.filterMask = other.filterMask
        self.outlierMap = other.outlierMap

DLSConfig.factory()

class MultiDataVector(DataVector):
    """Stores multiple values for each data point and can switch between a
    flat and a two-dimensional representation. For example, for multiple
    vectors consisting of N values, it stores *f(x)*, *g(x)* and *h(x)*
    (3 columns, N rows).
    The flat representation means to concatenate the columns of *f*, *g* and
    *h* (1 column, 3N rows) or repeating the first column *count* times.
    """
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
        rowCount = len(a) // self._count
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
                   "scatteringVector", "gammaDivR", "tauGamma",
                   "capTime" ,"countRate")

    @classproperty
    @classmethod
    def displayDataDescr(cls):
        return ("Title", "Measurements", "#Data",
                "Angle", "Description")

    @classproperty
    @classmethod
    def displayData(cls):
        return ("title", "measIndicesStr", "count",
                "anglesToStr", "description")

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

    def setMeasIndices(self, measIndices, a = None):
        """Sets the measurement index of this data set. Expects a list of
        tuples containing two integers of the measurement group index and the
        measurement index within that group.
        Stores the indices on a per angle basis which differs for each angle.
        """
        def verifyMeasIndex(measIndex):
            assert isList(measIndex) and len(measIndex) == 2
            assert all([isInteger(i) for i in measIndex])
        assert isList(measIndices)
        [verifyMeasIndex(mi) for mi in measIndices]
        # using tuples; setAngles() tests for tuple and hash(tuple()) works
        measIndices = tuple(measIndices)
        if isList(self.angles):
            try: # try to set the indices of a specific angle
                a = max(0, min(int(a), len(self.angles)-1))
                self._measIndices[a] = measIndices
                # sync length of measIndices with count of angles
                if len(self._measIndices) > len(self.angles):
                    del self._measIndices[len(self.angles):]
                elif len(self._measIndices) < len(self.angles):
                    self._measIndices[len(self._measIndices):len(self.angles)] = [()]
            except TypeError: # idx == None? set all angles to the same idx
                self._measIndices = [measIndices for i in range(len(self.angles))]
        else: # no angles yet
            self._measIndices = measIndices

    @property
    def measIndicesStr(self):
        # http://stackoverflow.com/questions/2154249
        def groupIndices(indices):
            """Identify groups of continuous numbers in a list."""
            ranges = []
            for key, group in groupby(enumerate(indices),
                                      lambda x: x[0] - x[1]):
                group = list(map(itemgetter(1), group))
                if len(group) > 1:
                    ranges.append("{0}-{1}".format(group[0], group[-1]))
                else:
                    ranges.append(str(group[0]))
            return ranges

        if self.measIndices is None:
            return ""
        measIndices = set(chain(*self.measIndices))
        summary = OrderedDict()
        for g, i in sorted(list(measIndices)):
            if g not in summary:
                summary[g] = []
            summary[g].append(i)
        # output format example: "0'1-3,5 1'2-4,6,8"
        # from group 0, indices 1,2,3,5 and from group 1, indices 2,3,4,6,8
        # https://msdn.microsoft.com/en-us/library/windows/desktop/aa365247.aspx
        res = " ".join(["{grp}'{indices}".format(grp = g,
            indices = ",".join(groupIndices(lst)))
                for g, lst in summary.items()])
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

    def _verifyDataMember(self, inArray, inArrayU = None):
        assert self.isValidInput(inArray), "Invalid data from file!"
        assert inArray.shape[1] == self.numAngles, (
            "Input data: {} columns != {} scattering angles"
            .format(inArray.shape[1], self.numAngles))
        if inArrayU is not None:
            assert self.isValidInput(inArrayU), "Invalid data from file!"
            assert inArray.shape == inArrayU.shape, \
                "Shapes of the provided data and its uncertainties do not match!"

    def setCorrelation(self, inArray, inArrayU = None):
        self._verifyDataMember(inArray, inArrayU)
        if inArrayU is None:
            inArrayU = zeros_like(inArray)
        self.f = MultiDataVector(u"G_2(τ)-1", inArray, rawU = inArrayU)

    def setCapTime(self, timeUnit, inArray):
        self._capTime = MultiDataVector("t", inArray.flatten(),
                                    unit = timeUnit, count = self.numAngles)

    def setCountRate(self, inArray, inArrayU = None):
        self._verifyDataMember(inArray, inArrayU)
        if inArrayU is None:
            inArrayU = zeros_like(inArray)
        self._countRate = MultiDataVector("Cnt(t)", inArray, rawU = inArrayU)

    @property
    def count(self):
        return len(self.x0.sanitized)

    # scattering angles

    def setAngles(self, anglesUnit, angles):
        """Expects angles in si units and their associated unit for proper
        formatting."""
        self._angles = angles.copy() # make them contiguous for hashing
        self._anglesUnit = anglesUnit
        # update measIndices using tuples, see setMeasIndices()
        # setAngles() tests for tuple and hash(tuple()) works
        measIndices = self.measIndices
        if isinstance(measIndices, tuple):
            # set the same measIndices for all angles
            self.setMeasIndices(measIndices)
        elif isinstance(measIndices, list):
            # contains already an idx tuple per angle
            for i in range(len(self.angles)):
                if i < len(measIndices):
                    self.setMeasIndices(measIndices[i], i)
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
            str(fmt).format(unit.toDisplay(a), unit.displayMagnitudeName)
                    for a in self.angles))

    @property
    def seriesKey(self):
        return "anglesToDisplay"

    @property
    def seriesKeyName(self):
        return (u"Scattering angle ({})"
                .format(self.anglesUnit.displayMagnitudeName))

    @property
    def anglesToDisplay(self):
        """Formatted scattering angles"""
        values = [self.anglesUnit.toDisplay(a) for a in self.angles]
        if len(values) == 1:
            values = values[0]
        return values

    # temperature, viscosity, refractiveIndex, wavelength including std.err

    def setTemperature(self, temp, stddev = 0.):
        self._temperature = (temp, stddev)
        self._calcGammaDivR()

    def setViscosity(self, vis, stddev = 0.):
        self._viscosity = (vis, stddev)
        self._calcGammaDivR()

    def setRefractiveIndex(self, refIdx, stddev = 0.):
        self._refractiveIndex = (refIdx, stddev)
        self._calcScatteringVector()

    def setWavelength(self, wavelen, stddev = 0.):
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
        self._gammaDivR = (  (- q2 * self.temperature[0] * KB)
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
        ddof = 1
        if len(others) == 1:
            ddof = 0 # prevent division by zero in numpy.std()
        # average basic properties
        for prop in ("temperature", "viscosity",
                     "refractiveIndex", "wavelength"):
            arr = array([getattr(o, prop)[0] for o in others])
            setFunc = getattr(self, _propSetterName(prop))
            if isCallable(setFunc):
                setFunc(arr.mean(), arr.std(ddof = ddof))
        # angles unchanged, but ensure they're identical everywhere
        assert array([self.angles == o.angles for o in others]).all(), \
               "Scattering angles differ between all DLS data to be combined!"
        # combine all measurement indices, usually one index for each o in others
        mis = set()
        for mi in [o.measIndices for o in others]:
            for i in chain(*mi):
                mis.add(i)
        self.setMeasIndices(tuple(sorted(list(mis))))
        # calculate average correlation values and their standard deviation
        stacked   = dstack([o.correlation.rawDataSrcShape for o in others])
        stackedCR = None
        if all([o.countRate is not None for o in others]):
            stackedCR = dstack([o.countRate.rawDataSrcShape for o in others])
        # filter outliers, possibly
        corr, corrU, cr, crU = None, None, None, None
        if (self.config.filterMask is not None
            and self.config.filterMask.shape == (len(self.angles), len(others))):
            # use prepared filter mask if shapes match
            corr  = np.zeros(stacked.shape[0:2])
            corrU = np.zeros(stacked.shape[0:2])
            if stackedCR is not None:
                cr    = np.zeros(stackedCR.shape[0:2])
                crU   = np.zeros(stackedCR.shape[0:2])
            for a, angle in enumerate(self.angles):
                # different count of good data for each angle average
                corrClean  = stacked[:,a,self.config.filterMask[a]]
                corr[:,a]  = corrClean.mean(-1)
                corrU[:,a] = corrClean.std(-1, ddof = ddof)
                if stackedCR is not None:
                    crClean  = stackedCR[:,a,self.config.filterMask[a]]
                    cr[:,a]  = crClean.mean(-1)
                    crU[:,a] = crClean.std(-1, ddof = ddof)
                if (self.config.outlierMap is None
                    or angle not in self.config.outlierMap
                    or not len(self.config.outlierMap[angle])):
                    continue
                mis = set(self.measIndices[a]) - self.config.outlierMap[angle]
                self.setMeasIndices(tuple(sorted(list(mis))), a)
                logging.warn("outliers for {a:.0f}{u}: {o}".format(
                                a = self.anglesUnit.toDisplay(angle),
                                u = self.anglesUnit.displayMagnitudeName,
                                o = sorted(list(self.config.outlierMap[angle]))))
        else:
            corr, corrU = stacked.mean(-1), stacked.std(-1, ddof = ddof)
            if stackedCR is not None:
                cr, crU     = stackedCR.mean(-1), stackedCR.std(-1, ddof = ddof)
        # combine the mean across all data sets with the existing tau
        # combine the std. deviation with the existing tau
        self.setCorrelation(corr, corrU)
        assert len(self.x0.siData) == len(self.f.siData), \
            "Dimensions of flattened data arrays do not match!"
        # same for count rate data
        if stackedCR is not None:
            self.setCountRate(cr, crU)
            assert len(self.capTime.siData) == len(self.countRate.siData), \
                "Dimensions of flattened data arrays do not match!"
        # reset config in order to fix callbacks, but maintain option values
        oldConfig = self.config
        self.initConfig()
        self.config.update(oldConfig)
        return self

    def splitPerAngle(self):
        lst = []
        if self.numAngles == 1:
            return [self]
        measIndices = copy.copy(self.measIndices)
        for i in range(self.numAngles):
            another = copy.copy(self)
            another.setAngles(self.anglesUnit, self.angles[i, newaxis])
            another.setMeasIndices(measIndices[i])
            another.setTau(self.tau.unit, self.tau.rawDataSrcShape)
            another.setCorrelation(
                    self.correlation.rawDataSrcShape[:, i, newaxis],
                    self.correlation.rawDataUSrcShape[:, i, newaxis])
            if isList(self.capTime):
                another.setCapTime(self.capTime.unit, self.capTime.rawDataSrcShape)
            if isList(self.countRate):
                another.setCountRate(self.countRate.rawDataSrcShape[:, i, newaxis],
                                     self.countRate.rawDataUSrcShape[:, i, newaxis])
            # reset config in order to fix callbacks, but maintain option values
            another.initConfig()
            another.config.update(self.config) # forward settings
            lst.append(another)
        return lst

    @staticmethod
    def preProcess(dataList): # FIXME: rename, e.g. to 'averageByAngle'
        """Starts accumulation of related data sets among the currently loaded
        ones. Finally, removes such source data sets and adds the new combined
        one."""
        if not isIterable(dataList):
            return # nothing to do
        samples = OrderedDict()
        def makeKey(data):
            key = (data.title,) # == self.sampleName in CGSFile
            if hasattr(data, "angles"):
                key += tuple(data.angles)
            return key
        for d in dataList: # group data objects by their title
            key = makeKey(d)
            if key not in samples:
                samples[key] = []
            samples[key].append(d)
        if not len(samples):
            return # nothing to do
        dataList = [] # accumulate data objects if possible
        for dummy, lst in samples.items():
            if not len(lst):
                continue
            if hasattr(lst[0], "analyseOutliers"):
                lst[0].analyseOutliers(lst)
            avg = None
            if (hasattr(lst[0], "accumulate")
                and hasattr(lst[0].config, "doAverage")
                and lst[0].config.doAverage()):
                avg = lst[0].accumulate(lst)
            if avg is None:
                dataList.extend(lst)
            else:
                dataList.append(avg)
        # add the combined dls data split up per angle
        def splitUp(d):
            return d.splitPerAngle()
            #try: # perhaps test for isinstance(d, DLSData) instead
            #    return d.splitPerAngle()
            #except AttributeError:
            #    return (d,)
        dataList = [s for dl in (splitUp(d) for d in dataList) for s in dl]
        return dataList

    def analyseOutliers(self, dataList):
        if dataList is None or len(dataList) < 2:
            return
        if any([d.countRate is None for d in dataList]):
            return # None or inconsistent data to analyse
        # outlier detection by median absolute deviation from:
        # http://stackoverflow.com/a/22357811
        stacked = dstack((d.countRate.rawDataSrcShape for d in dataList))
        # median over stacking axis, extend array shape by one dim
        median = np.median(stacked, axis = -1)[:,:,None]
        diff = np.sqrt(np.sum((stacked - median)**2, axis = 0)) # sum along tau
        self.config.filterMask = np.ones(diff.shape, dtype = bool)
        self.config.medianCountRate = dict()
        self.config.outlierMap = dict()
        # median along observations, for each angle
        medAbsDev = np.median(diff, axis = -1)[:,None]
        if any(medAbsDev == 0.):
            return
        modifiedZScore = 0.6745 * diff / medAbsDev
        goodData = modifiedZScore < self.config.outlierThreshold()
        # store the median for plotting later
        self.config.filterMask = goodData
        for a, angle in enumerate(self.angles):
            self.config.medianCountRate[angle] = np.squeeze(median[:,a])
            # for each angle, store illegal measIndices
            self.config.outlierMap[angle] = set([d.measIndices[a][0]
                for i, d in enumerate(dataList) if not goodData[a][i]])

    def ids(self):
        """Returns a text containing the objects ids of the embedded data
        objects. (for debugging purposes)"""
        out = [repr(self)]
        for p in self._properties:
            out.append("{}: {}".format(p, id(getattr(self, _privPropName(p)))))
        return u"\n".join(out)

    def __init__(self, **kwargs):
        super(DLSData, self).__init__(**kwargs)

    def setConfig(self, config = None):
        """Overrides DataObj.setConfig() to forward some general settings to
        all configurations before the sampleName filter kicks in."""
        if not isinstance(config, self.configType):
            return # ignore configurations of other types
        if self.config is not None:
            self.config.doAverage.setValue(config.doAverage())
            self.config.plotCountRate.setValue(config.plotCountRate())
        super(DLSData, self).setConfig(config)

    def updateConfig(self):
        # this has to be set before base class updateConfig()
        self.config.x0Low.setUnit(self.tau.unit)
        self.config.x0High.setUnit(self.tau.unit)
        super(DLSData, self).updateConfig()

    def __str__(self):
        out = [u"## {0} '{1}'".format(self.__class__.__name__, self.title)]
        for p in self._properties:
            if "anglesUnit" in p: continue
            value = getattr(self, _privPropName(p))
            if "angles" in p: value = self.anglesToStr
            out.append(u"{0}: {1}".format(p, value))
        return u"\n".join(out)

    def __hash__(self):
        value = hash(self.title) ^ hash(self.filename)
        for p in self._properties:
            propData = getattr(self, _privPropName(p))
            if isinstance(propData, list): # fixes measIndices hashing
                propData = tuple(propData)
            try:
                value ^= hash(propData)
            except TypeError: # numpy.ndarray
                value ^= hashNumpyArray(propData)
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
