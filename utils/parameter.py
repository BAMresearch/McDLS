# -*- coding: utf-8 -*-
# utils/parameter.py

from utils import mixedmethod, isList, testfor, isInteger
from cutesnake.algorithm import (ParameterBase, ParameterFloat,
                                 ParameterNumerical, ParameterBoolean,
                                 ParameterLog, ParameterString)
from cutesnake.algorithm import Parameter
from cutesnake.dataset import DataSet, DisplayMixin

import logging
import numpy as np
from cutesnake.utils.tests import testfor
from cutesnake.algorithm.numbergenerator import NumberGenerator, RandomUniform  
from cutesnake.algorithm.parameter import ParameterError, ValueRangeError

class ParameterGeneratorError(ParameterError):
    pass

def _makeProperty(varName):
    def getter(selforcls):
        return getattr(selforcls, varName)
    return property(getter)

class Moments(object):
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

    def __init__(self, contribs, paramIndex, valueRange, fraction, algo = None):
        self._setValidRange(contribs[:, paramIndex, :], valueRange)
        self._calcMoments(contribs[:, paramIndex, :], fraction)
        if algo is not None:
            scalingFactors = algo.result[paramIndex]['scalingFactors']
            self._calcPartialIntensities(contribs, scalingFactors, algo)

    def _setValidRange(self, contribs, valueRange):
        """Calculate contributions mask to be within the given value range.
        contribs: parameter value sets of shape
                  <number of contributions x number of repetitions>
        """
        testfor(contribs.ndim == 2, ValueError)
        numContribs, numReps = contribs.shape
        self._validRange = np.zeros_like(contribs.T, dtype = bool)
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
        val = np.zeros(numReps)
        mu  = np.zeros(numReps)
        var = np.zeros(numReps)
        skw = np.zeros(numReps)
        krt = np.zeros(numReps)
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
            sigma   = np.sqrt(abs(var[ri]))
            skw[ri] = ( sum( (rset-mu[ri])**3 * frac )
                     / (sum(frac) * sigma**3))
            krt[ri] = ( sum( (rset-mu[ri])**4 * frac )
                     / (sum(frac) * sigma**4))

        self._total    = (val.mean(), val.std(ddof = 1))
        self._mean     = ( mu.mean(),  mu.std(ddof = 1))
        self._variance = (var.mean(), var.std(ddof = 1))
        self._skew     = (skw.mean(), skw.std(ddof = 1))
        self._kurtosis = (krt.mean(), krt.std(ddof = 1))

    #partial intensities not parameter-specific. Move to model/histogram?
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
        partialIntensities = np.zeros((numReps, data.q.shape[0]))
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

class VectorResult(object):
    """Stores multiple populations of a single result data vector.
    Calculates statistics at initialization."""
    _mean = None
    _std = None
    _full = None

    # some read-only attributes

    @property
    def mean(self):
        return self._mean

    @property
    def std(self):
        return self._std

    @property
    def full(self):
        return self._full

    def __init__(self, vecResult):
        assert vecResult.ndim == 2 # 2 dim input
        self._full = vecResult
        self._mean = self._full.mean(axis = 1)
        self._std = self._full.std(axis = 1)

# put this sketch here, for the moment can be placed in a separate file later
class Histogram(DataSet, DisplayMixin):
    """Stores histogram related settings of a parameter.
    The results too, eventually(?). yes, please.
    Stores&calculates rangeInfo() results for all available weighting options.
    """
    _param      = None # the FitParameter this histogram belongs to
    _binCount   = None # list of bin counts
    _xscale     = None # list of scalings
    _yweight    = None # list of weightings
    _xrange     = None # list of tuples/pairs
    _stats      = None # rangeInfo() results, RangeStats lists

    # results: class HistogramResult?
    _xLowerEdge = None
    _xMean = None
    _xWidth = None
    _bins = None
    _cdf = None
    _observability = None
    _moments = None

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

    @property
    def xLowerEdge(self):
        return self._xLowerEdge

    def _setXLowerEdge(self):
        # Now bin whilst keeping track of which contribution ends up in
        # which bin: set bin edge locations
        if 'lin' in self.xscale:
            # histogramXLowerEdge contains #histogramBins+1 bin edges,
            # or class limits.
            self._xLowerEdge = np.linspace(
                    self.lower, self.upper, self.binCount + 1)
        else:
            self._xLowerEdge = np.logspace(
                    np.log10(self.lower), np.log10(self.upper),
                    self.binCount + 1)
        self._setXWidth()
        self._setXMean()

    @property
    def xMean(self):
        return self._xMean

    def _setXMean(self):
        if self.xLowerEdge is None:
            return
        # NOTE: isn't this the same as:
        # self.xWidth * .5 + self.xLowerEdge
        self._xMean = np.zeros(self.binCount)
        for i in range(self.binCount):
            self._xMean[i] = self.xLowerEdge[i:i+2].mean()

    @property
    def xWidth(self):
        return self._xWidth

    def _setXWidth(self):
        if self.xLowerEdge is None:
            return
        self._xWidth = np.diff(self.xLowerEdge)

    @property
    def bins(self):
        return self._bins

    @property
    def cdf(self):
        return self._cdf

    @property
    def observability(self):
        return self._observability

    def _setObservability(self, allObservability):
        self._observability = np.zeros(self.binCount)
        if allObservability is None:
            return
        testfor(allObservability.shape[0] == self.binCount, ValueError)
        for bi in range(self.binCount):
            # for observabilities over all repetitions select the largest
            obs = allObservability[bi, :]
            self._observability[bi] = obs[obs < np.inf].max()

    @property
    def moments(self):
        return self._moments

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

    def _binMask(self, bini, parValues):
        # indexing which contributions fall into the radius bin
        return (  (parValues >= self.xLowerEdge[bini])
                * (parValues <  self.xLowerEdge[bini + 1]))

    def calc(self, contribs, paramIndex, fractions):
        self._setXLowerEdge()
        self._calcRepetitions(contribs, paramIndex, fractions)

    def _calcRepetitions(self, contribs, paramIndex, fractions):
        numContribs, dummy, numReps = contribs.shape
        binLst, obsLst, cdfLst = [], [], []
        fractions, minReq = fractions[self.yweight]
        for ri in range(numReps):
            parValues = contribs[:, paramIndex, ri]
            bins, binObs, cdf = self._calcBins(
                    contribs, parValues, fractions[:, ri], minReq[:, ri])
            binLst.append(bins)
            obsLst.append(binObs)
            cdfLst.append(cdf)
        # set final result: y values, CDF and observability of all bins
        self._bins = VectorResult(np.vstack(binLst).T)
        self._cdf = VectorResult(np.vstack(cdfLst).T)
        self._setObservability(np.vstack(obsLst).T)
        self._moments = Moments(contribs, paramIndex, self.xrange, fractions)

    def _calcBins(self, contribs, parValues, fraction, minReq):
        """Returns np arrays for the bin values, observability and the CDF
        based on the bin values."""
        # single set of R for this calculation
        bins = np.zeros(self.binCount)
        binObs = np.zeros(self.binCount)
        for bi in range(self.binCount):
            val, obs = self._calcBin(
                    self._binMask(bi, parValues),
                    fraction, minReq)
            bins[bi] = val
            binObs[bi] = obs
            cdf = self._calcCDF(bins)
        return bins, binObs, cdf

    def _calcBin(self, binMask, fraction, minReq):
        """
           *fraction*: overall fraction (number or volume, weighting depending)
        """
        # y contains the volume fraction for that radius bin
        binValue = sum(fraction[binMask])
        if np.isnan(binValue):
            binValue = 0.
        # observability below
        if not any(binMask):
            binMinReq = 0.
        else:
            binMinReq = minReq[binMask].mean()
        return binValue, binMinReq

    def _calcCDF(self, bins):
        cdf = np.zeros_like(bins)
        cdf[0] = bins[0]
        for i in range(1, len(cdf)):
            cdf[i] = cdf[i - 1] + bins[i]
        return cdf

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
        out = ["hist"]
        for attr in self.displayData:
            val = getattr(self, attr)
            out.append(str(val))
        idx = self.displayData.index("paramName")
        # replace the parameter display name by its internal name
        # (short, no spaces)
        out[idx+1] = self.param.name()
        return "-".join(out)

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
        logging.debug('Hist init: {} [{}, {}]'
                .format(param.name(), lower, upper))
        DataSet.__init__(self, "({0}, {1})".format(lower, upper), None)
        # add it to the parameter here
        self.param = param # parameter we belong to is mandatory
        self.binCount = binCount # bin count is mandatory
        self.xrange = (lower, upper)
        # setter chose the first option available for invalid options
        self.xscale = xscale
        self.yweight = yweight

    # TODO: Function toString() or toJson() or similar which makes it
    # serializable and also a classmethod which constructs it from serial data

class Histograms(list):
    """Manages a set of user configured histograms for evaluation after
    monte-carlo run."""
    def append(self, value):
        if value in self:
            return
        list.append(self, value)

    def updateRanges(self):
        # work directly on the histograms, no copy
        for i in range(len(self)):
            self[i].updateRange()

    def calc(self, *args):
        # work directly on the histograms, no copy
        for i in range(len(self)):
            self[i].calc(*args)

def isActiveParam(param):
    """Checks any type of parameter for activeness.
    Shorter than that below or a try/except clause."""
    return isinstance(param, FitParameterBase) and param.isActive()

class FitParameterBase(ParameterBase):
    """Deriving parameters for curve fitting from
    cutesnake.algorithm.parameter to introduce more specific fit
    related attributes."""
    ParameterBase.setAttributes(locals(), histograms = None,
            activeValues = list(), activeRange = (None, None))

    def __init__(self):
        super(FitParameterBase, self).__init__()
        # point the parameter reference to this instance now
        # (instead of the type previously)
        if self.isActive():
            for i in range(len(self.histograms())):
                self.histograms()[i].param = self

    @mixedmethod
    def setValueRange(selforcls, newRange):
        try:
            super(FitParameterBase, selforcls).setValueRange(newRange)
            selforcls.histograms().updateRanges()
        except:
            pass

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
            if None in selforcls.activeRange():
                lo, hi = selforcls.valueRange()
            else:
                lo, hi = selforcls.activeRange()
            selforcls.setHistograms(Histograms()) # init histogram list
            # add a default histogram
            # NOTE: apply default bin count from MCSASParameters here
            # or in Histogram.__init__ defaults
            selforcls.histograms().append(Histogram(selforcls, lo, hi))
        elif not isActive:
            selforcls.setHistograms(None)

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
    
    # Following should be moved to FitParameterNumerical
    ParameterBase.setAttributes(locals(), "generator")

    @mixedmethod
    def setActiveRange(selforcls, newRange):                                       
        # tests nicked from above 
        testfor(len(newRange) == 2, ValueRangeError,
                "Active ranges have to consist of two values!") 
        # always clip range settings to min/max values:
        newRange = selforcls.clip(newRange)
        # sets range for active fitting parameter limits
        selforcls._activeRange = (min(newRange), max(newRange))

    @mixedmethod
    def activeRange(selforcls): 
        if selforcls._activeRange is None:
            return (None, None)
        else:
            return selforcls._activeRange

    @mixedmethod                                                               
    def displayActiveRange(selforcls):                                         
        """value bounds in display units used for parameter generator"""       
        magConv = selforcls.unit.magnitudeConversion()                         
        vRange = selforcls.activeRange()                                       
        newRange = (min(vRange) / magConv, max(vRange) / magConv)              
        return newRange                                                        

    @mixedmethod                                                               
    def setDisplayActiveRange(selforcls, newRange):                            
        """sets value range after converting input from display to SI units""" 
        magConv = selforcls.unit.magnitudeConversion()                         
        newRange = (min(newRange) * magConv, max(newRange) * magConv)          
        selforcls.setActiveRange(newRange)                                     

    @mixedmethod                                                               
    def setGenerator(selforcls, newGenerator):
        if isinstance(newGenerator, type):
            testfor(issubclass(newGenerator, NumberGenerator), 
                    ParameterGeneratorError, "NumberGenerator type expected!")
        else:
            newGenerator = RandomUniform
        selforcls._generator = newGenerator

    def generateValues(selforcls, numberGenerator, 
            defaultRange, lower, upper, count):        
        # works with vectors of multiple bounds too                                
        vRange = defaultRange
        if lower is None:
            lower = vRange[0]
        if upper is None:
            upper = vRange[1] 
        vRange = (np.maximum(vRange[0], lower), np.minimum(vRange[1], upper))
        if isList(vRange[0]) and isList(vRange[1]):
            assert len(vRange[0]) == len(vRange[1]), "Provided value range is unsymmetrical!"
        try: # update count to length of provided bound vectors                    
            count = max(count, min([len(x) for x in vRange]))                      
        except:                                                                    
            pass                                                                   
        values = numberGenerator.get(count)                                        
        # scale numbers to requested range                                         
        return values * (vRange[1] - vRange[0]) + vRange[0]                        

    def generate(self, lower = None, upper = None, count = 1):                 
        return self.generateValues(self.generator(), self.activeRange(), 
                lower, upper, count).astype(self.dtype)  

class FitParameterString(FitParameterBase, ParameterString):
    pass

class FitParameterBoolean(FitParameterBase, ParameterBoolean):
    pass

class FitParameterNumerical(FitParameterBase, ParameterNumerical):

    def __str__(self):
        lo, hi = self.displayActiveRange()
        displayValue = u', active range: [{0}, {1}]'.format(lo, hi)
        return (super(FitParameterNumerical, self).__str__() + u"{0} ({1})".format(
                displayValue,
                self.suffix()))

class FitParameterFloat(FitParameterNumerical, ParameterFloat):
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
