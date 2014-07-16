# -*- coding: utf-8 -*-
# utils/parameter.py

from utils import mixedmethod, isList, testfor, isInteger
from cutesnake.algorithm import (ParameterBase, ParameterFloat,
                                 ParameterNumerical, ParameterBoolean,
                                 ParameterLog, ParameterString)
from cutesnake.algorithm import Parameter
from cutesnake.dataset import DataSet, DisplayMixin

import logging
import numpy
from cutesnake.utils.tests import testfor
def _makeProperty(varName):
    def getter(selforcls):
        return getattr(selforcls, varName)
    return property(getter)

class RangeStats(object):
    _intensity  = None # partial intensity
    _total      = None
    _mean       = None # 1st moment
    _variance   = None # 2nd moment
    _skew       = None # 3rd moment
    _kurtosis   = None # 4th moment
    # intermediate valid contribs indices for the given parameter range
    _validRange = None

    # create getter properties for scalar attributes
    for varName in ("_intensity", "_total", "_mean", "_variance", "_skew",
                    "_kurtosis"):
        locals()[varName[1:]] = _makeProperty(varName)

    @staticmethod
    def fieldNames():
        """Returns the field names in the same order as str()"""
        return (
                "totalValue", "totalValueStd",
                      "mean", "meanStd",
                  "variance", "varianceStd",
                      "skew", "skewStd",
                  "kurtosis", "kurtosisStd"
                 )

    @property
    def fields(self):
        """Tuple of member data incl. uncertainty for export."""
        return (self._total + self._mean + self._variance
                + self._skew + self._kurtosis)

    def __str__(self):
        return str(self.fields)

    def __repr__(self):
        return str(self)

    def __init__(self, paramIndex, valueRange, fraction, algo):
        contribs = algo.result[0]['contribs']
        self._setValidRange(contribs[:, paramIndex, :], valueRange)
        self._calcMoments(contribs[:, paramIndex, :], fraction)
        #for par in algo.activeParams():
        #self._calcMoments(self.activeVal(), fraction)
        scalingFactors = algo.result[paramIndex]['scalingFactors']
        self._calcPartialIntensities(contribs, scalingFactors, algo)

    def _setValidRange(self, contribs, valueRange):
        """Calculate contributions mask to be within the given value range.
        contribs: parameter value sets of shape
                  <number of contributions x number of repetitions>
        """
        testfor(contribs.ndim == 2, ValueError)
        numContribs, numReps = contribs.shape
        self._validRange = numpy.zeros_like(contribs.T, dtype = bool)
        for ri in range(numReps):
            # the single set of R for this calculation
            rset = contribs[:, ri]
            self._validRange[ri] = ((rset > min(valueRange))
                                  * (rset < max(valueRange)))

    def _calcMoments(self, contribs, fraction):
        """Calculates the moments of the distribution of the current
        particular (implied) parameter.
        contribs: parameter value sets of shape
                  <number of contributions x number of repetitions>
        fraction: number or volume fraction
        """
        numContribs, numReps = contribs.shape
        val = numpy.zeros(numReps)
        mu  = numpy.zeros(numReps)
        var = numpy.zeros(numReps)
        skw = numpy.zeros(numReps)
        krt = numpy.zeros(numReps)
        # loop over each repetition
        for ri in range(numReps):
            # the single set of R for this calculation
            if not any(self._validRange[ri]):
                continue # what to do if validRange is empty?
            rset = contribs[self._validRange[ri], ri]
            frac = fraction[self._validRange[ri], ri]
            val[ri] = sum(frac)
            mu[ri]  = sum(rset * frac)/sum(frac)
            var[ri] = sum( (rset-mu[ri])**2 * frac )/sum(frac)
            sigma   = numpy.sqrt(abs(var[ri]))
            skw[ri] = ( sum( (rset-mu[ri])**3 * frac )
                     / (sum(frac) * sigma**3))
            krt[ri] = ( sum( (rset-mu[ri])**4 * frac )
                     / (sum(frac) * sigma**4))

        self._total    = (val.mean(), val.std(ddof = 1))
        self._mean     = ( mu.mean(),  mu.std(ddof = 1))
        self._variance = (var.mean(), var.std(ddof = 1))
        self._skew     = (skw.mean(), skw.std(ddof = 1))
        self._kurtosis = (krt.mean(), krt.std(ddof = 1))

    #partial intensities not parameter-specific. Move to model?
    def _calcPartialIntensities(self, contribs, scalingFactors, algo):
        """scalingFactors: scaling and background for each repetition.
        contribs: all contributions also of other parameters; not required
        once the results/contribs are stored in the model/parameters then
        algo.calcModel() does not need them provided externally
        """
        # now we can calculate the intensity contribution by the subset of
        # spheres highlighted by the range:
        numContribs, numParams, numReps = contribs.shape
        data = algo.dataPrepared
        # loop over each repetition
        partialIntensities = numpy.zeros((numReps, data.q.shape[0]))
        # Intensity scaling factors for matching to the experimental
        # scattering pattern (Amplitude A and flat background term b,
        # defined in the paper)
        for ri in range(numReps):
            if not any(self._validRange[ri]):
                continue # what to do if validRange is empty?
            rset = contribs[self._validRange[ri], :, ri]
            # compensated volume for each sphere in the set
            it, vset = algo.calcModel(data, rset)
            it = algo.model.smear(it)
            # a set of volume fractions
            partialIntensities[ri, :] = it.flatten() * scalingFactors[0, ri]

        self._intensity = (partialIntensities.mean(axis = 0),
                           partialIntensities.std(axis = 0))

# put this sketch here, for the moment can be placed in a separate file later
class Histogram(DataSet, DisplayMixin):
    """Stores histogram related settings of a parameter.
    The results too, eventually(?). yes, please.
    Stores&calculates rangeInfo() results for all available weighting options.
    """
    # back reference of the FitParameter this histogram belongs to
    _param      = None # this is not necessary here. belongs to one par.
    _binCount   = None # list of bin counts
    _xscale     = None # list of scalings
    _yweight    = None # list of weightings
    _xrange     = None # list of tuples/pairs
    _stats      = None # rangeInfo() results, RangeStats lists

    @property
    def param(self):
        return self._param

    @param.setter
    def param(self, newParam):
        self._param = newParam

    @property
    def paramName(self):
        return self.param.displayName()

    @property
    def binCount(self):
        return self._binCount

    @binCount.setter
    def binCount(self, binCount):
        testfor(isInteger(binCount) and binCount > 0,
                ValueError,
                "Histogram bin count has to be an integer larger 0!")
        self._binCount = max(0, int(binCount))

    @property
    def xscale(self):
        return self._xscale

    @xscale.setter
    def xscale(self, kind):
        """Sets it to the first available option by default."""
        self._xscale = str(kind).strip() # remove whitespace eventually
        if self._xscale not in self.xscaling():
            self._xscale = self.xscaling(0)

    @property
    def yweight(self):
        return self._yweight

    @yweight.setter
    def yweight(self, kind):
        self._yweight = str(kind).strip() # remove whitespace eventually
        if self._yweight not in self.yweighting():
            self._yweight = self.yweighting(0)

    @property
    def xrange(self):
        return self._xrange

    @xrange.setter
    def xrange(self, valueRange):
        lo, hi = min(valueRange), max(valueRange)
        # restrict the provided range to the current range of the parameter
        lo, hi = max(self.param.min(), lo), min(self.param.max(), hi)
        self._xrange = lo, hi

    @property
    def lower(self):
        return self._xrange[0]

    @property
    def upper(self):
        return self._xrange[1]

    def updateRange(self):
        """Updates histogram range according to a changed parameter range
        to keep it valid."""
        self.xrange = self.xrange # call getter & setter again

    @staticmethod
    def displayDataDescr():
        return (
                "parameter",
                "lower",
                "upper",
                "Number of bins",
                "X-axis scaling",
                "Y-axis weighting"
                )

    @property
    def displayData(self):
        return (
                "paramName",
                "lower",
                "upper",
                "binCount",
                "xscale",
                "yweight"
                )

    @staticmethod
    def xscaling(index = None):
        avail = ('lin', 'log')
        try:
            return avail[index]
        except:
            return avail

    @staticmethod
    def yweighting(index = None):
        avail = ('vol', 'num')
        try:
            return avail[index]
        except:
            return avail

    def addRange(self, lower, upper):
        """Adds the provided range and returns its index.
        The index may be different from the last if the range exists."""
        newRange = (min(lower, upper), max(lower, upper))
        try:
            return self.ranges.index(newRange)
        except ValueError:
            self._ranges.append(newRange)
        return len(self._ranges) - 1

    def resetRanges(self):
        self._ranges = []

    def calcStats(self, paramIndex, algo):
        """Calculates distribution statistics for all available weighting
        options over all ranges previously set."""
        self._stats = [] # statistics for all ranges
        for valueRange in self.ranges:
            allWeighting = [] # for all weighting options of a range
            for weighting in self.weighting():
                if weighting == 'vol':
                    fraction = algo.result[paramIndex]['volumeFraction']
                elif weighting == 'num':
                    fraction = algo.result[paramIndex]['numberFraction']
                logging.info("Calculating {weighting} weighted distribution "
                             "statistics of {param} within {range} ..."
                             .format(weighting = weighting.upper(),
                                     param = self.param.name(),
                                     range = valueRange))
                allWeighting.append(RangeStats(paramIndex, valueRange,
                                                fraction, algo))
            self._stats.append(allWeighting)
        return self._stats

    @property
    def stats(self):
        return self._stats

    def iterStats(self):
        for valueRange, rangeStats in zip(self.ranges, self.stats):
            for weighting, weightingStats in zip(self.weighting(), rangeStats):
                yield valueRange, weighting, weightingStats

    def __str__(self):
        out = []
        for attr in self.displayData:
            val = getattr(self, attr)
            out.append(str(val))
        return "Hist(" + ", ".join(out) + ")"

    __repr__ = __str__

    def __eq__(self, other):
        """Compares with another Histogram or tuple."""
        if id(self.param) != id(other.param):
            return False
        for attr in "lower", "upper", "binCount", "xscale", "yweight":
            thisval = getattr(self, attr)
            otherval = getattr(other, attr)
            if otherval != thisval:
                return False
        return True

    def __neq__(self, other):
        return not self.__eq__(other)

    def __init__(self, param, lower, upper, binCount = 50,
                 xscale = None, yweight = None):
        """Creates an histogram with default bin count, will be updated later."""
        DataSet.__init__(self, "({0}, {1})".format(lower, upper), None)
        # add it to the parameter here
        self.param = param # parameter we belong to is mandatory
        self.binCount = binCount # bin count is mandatory
        self.xrange = (lower, upper)
        # setter chose the first option available for invalid options
        self.xscale = None
        self.yweight = None
        print " Histogram:", lower, upper, self.xrange, self.xscale, self.yweight

    # TODO: Function toString() or toJson() or similar which makes it
    # serializable and also a classmethod which constructs it from serial data

class Histograms(list):
    """Manages a set of user configured histograms for evaluation after
    monte-carlo run."""
    def append(self, value):
        list.append(self, value)

    def updateRanges(self):
        for i in range(0, len(self)):
            self[i].updateRange()

class FitParameterBase(ParameterBase):
    """Deriving parameters for curve fitting from
    cutesnake.algorithm.parameter to introduce more specific fit
    related attributes."""
    # by default it is not fitted, inactive
    ParameterBase.setAttributes(locals(), histograms = None,
                                activeValues = list())

    def __init__(self):
        super(FitParameterBase, self).__init__()
        # point the parameter reference to this instance now
        # (instead of the type previously)
        if self.isActive():
            for i in range(0, len(self.histograms())):
                self.histograms()[i].param = self

    @mixedmethod
    def isActive(selforcls):
        """Tests if there is an histogram defined.
        Provided for convenience."""
        return isList(selforcls.histograms())

    @mixedmethod
    def setActive(selforcls, isActive):
        """Sets this parameter as active. It will be fitted.
        Temporary in replacement of setIsActive(). Can be removed later if
        histgram setup for each parameter is implemented elsewhere"""
        if isActive and not selforcls.isActive():
            # set only if there is no histogram defined already
            print "setactive", selforcls
            lo, hi = selforcls.valueRange()
            selforcls.setHistograms(Histograms()) # init histogram list
            # add a default histogram
            # NOTE: apply default bin count from MCSASParameters here
            # or in Histogram.__init__ defaults
            selforcls.histograms().append(Histogram(selforcls, lo, hi))
        elif not isActive:
            selforcls.setHistograms(None)

    @mixedmethod
    def setValueRange(selforcls, newRange):
        try:
            super(FitParameterBase, selforcls).setValueRange(newRange)
            selforcls.histograms().updateRanges()
        except:
            pass

    @mixedmethod
    def activeVal(selforcls, val, index = None):
        """
        activeVal is set after a successful MC run. It is a list of arrays, 
        with each array storing the parameter values of a successful run.
        If index is supplied, only the array at that list index is returned, 
        otherwise the entire list is returned.
        """
        if index is None:
            return selforcls.activeValues()
        else:
            return selforcls.activeValues()[
                    index%len(selforcls.activeValues())]

    @mixedmethod
    def setActiveVal(selforcls, val, index = None):
        """
        Sets or appends to the list of activeVal. There *could* be a race
        condition if two mcFit instances try to append at the same time. 
        Therefore, an index can be supplied to identify the list index to set 
        or change. Each mcFit instance should be assigned a serial (index) 
        number.
        If the list is not long enough to accommodate the value, it will be 
        extended.
        """
        if index is None:
            # append to end
            index = len(selforcls.activeValues())

        while len(selforcls.activeValues()) < (index + 1):
            # expand list to allow storage of value
            tempVal = selforcls.activeValues()
            tempVal.append(None)
            selforcls.setActiveValues(tempVal)

        if not selforcls.isActive():
            logging.error(
            'value of parameter cannot be set, parameter not active')
            return

        tempVal = selforcls.activeValues()
        tempVal[index] = val
        selforcls.setActiveValues(tempVal)

class FitParameterString(FitParameterBase, ParameterString):
    pass

class FitParameterBoolean(FitParameterBase, ParameterBoolean):
    pass

class FitParameterNumerical(FitParameterBase, ParameterNumerical):
    pass

class FitParameterFloat(FitParameterBase, ParameterFloat):
    pass

class FitParameterLog(FitParameterBase, ParameterLog):
    pass

def FitParameter(*args, **kwargs):
    return Parameter(*args, paramTypes = (
                FitParameterBoolean,
                FitParameterFloat,
                FitParameterNumerical,
                FitParameterString,
                FitParameterBase
            ), **kwargs)

# vim: set ts=4 sts=4 sw=4 tw=0:
