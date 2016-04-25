# -*- coding: utf-8 -*-
# scattering/dataobj.py

"""
Represents input data associated with a measurement.
"""

from __future__ import absolute_import # PEP328
import os # Miscellaneous operating system interfaces
from numpy import all as np_all
import numpy as np

from abc import ABCMeta, abstractproperty, abstractmethod
from bases.dataset import DataSet, DisplayMixin
from utils.units import Unit, NoUnit
from utils import classproperty
import logging

class DataVector(object):
    """ a class for combining aspects of a particular vector of data.
    This is intended only as a storage container without additional functionality.
    """
    _name = None # descriptor for axes and labels, unicode string at init.
    _raw = None # Relevant data directly from input file, non-SI units, unsanitized
    _rawU = None # raw uncertainties, unsanitized
    _siData = None # copy raw data in si units, often accessed
    _siDataU = None # uncertainties in si units, sanitized, with eMin taken into account
    _unit = None # instance of unit
    _limit = None # two-element vector with min-max
    _validIndices = None # valid indices. 
    
    def __init__(self, name, raw, rawU = None, unit = None):
        self._name = name
        self._raw = raw
        self._rawU = rawU
        self.unit = unit
        self.validIndices = np.arange(self.raw.size) # sets limits as well

    @property
    def name(self):
        return unicode(self._name)

    @property
    def validIndices(self):
        return self._validIndices

    @validIndices.setter
    def validIndices(self, indices):
        assert indices.min() >= 0
        assert indices.max() <= self.siData.size
        self._validIndices = indices
        self._limit = [self.sanitized.min(), self.sanitized.max()]

    @property
    def sanitized(self):
        return self.siData.copy()[self.validIndices]

    @sanitized.setter
    def sanitized(self, val):
        assert(val.size == self.validIndices.size)
        self.siData[self.validIndices] = val

    @property
    def sanitizedU(self):
        return self.siDataU.copy()[self.validIndices]

    @sanitizedU.setter
    def sanitizedU(self, val):
        assert(val.size == self.validIndices.size)
        self.siDataU[self.validIndices] = val

    # siData
    @property
    def siData(self):
        return self._siData

    @siData.setter
    def siData(self, vec):
        self._siData = vec

    # siDataU, uncertainties on siData
    @property
    def siDataU(self):
        return self._siDataU
    
    @siDataU.setter
    def siDataU(self, vec):
        self._siDataU = vec

    # raw data values
    @property
    def raw(self):
        return self._raw

    # raw uncertainties
    @property
    def rawU(self):
        return self._rawU

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, newUnit):
        if not isinstance(newUnit, Unit):
            self._unit = NoUnit()
            self.siData = self.raw.copy()
            if self.rawU is not None:
                self.siDataU = self.rawU.copy()
        else:
            self._unit = newUnit
            self.siData = self.unit.toSi(self.raw)
            if self.rawU is not None:
                self.siDataU = self.unit.toSi(self.rawU)

    # TODO: define min/max properties for convenience?
    @property
    def limit(self):
        return self._limit

    @property
    def limsString(self):
        return u"{0:.3g} ≤ {valName} ({magnitudeName}) ≤ {1:.3g}".format(
                self.unit.toDisplay(self.limit[0]),
                self.unit.toDisplay(self.limit[1]),
                magnitudeName = self.unit.displayMagnitudeName,
                valName = self.name)

class DataObj(DataSet, DisplayMixin):
    """General container for data loaded from file. It offers specialised
    methods to derive information from the provided data.
    """
    __metaclass__ = ABCMeta
    _filename = None
    _config = None
    _validIndices = None # based on masks to filter certain values
    _validMask = None
    _x0 = None
    _x1 = None
    _x2 = None
    _f  = None
    _fBin = None # let's try to do it this way. 
    _x0Bin = None

    def reBin(self):
        """ 
        rebinning method, to be run (f.ex.) upon every "Start" buttonpress. 
        For now, this will rebin using the x0 vector as a base, although the 
        binning vector can theoretically be chosen freely.
        """
        logging.info("Initiating binning procedure")
        nBin = self.config.nBin.value()
        # self._binned = DataVector() once binning finishes.. dataVector can be set once.
        sanX = self.x0.sanitized
        x0Bin = np.zeros(nBin)
        fBin  = np.zeros(nBin)
        fuBin  = np.zeros(nBin)
        validMask = np.zeros(nBin, dtype = bool) #default false

        if not(nBin > 0):
            return # no need to do the actual rebinning. values stay None.

        # prepare bin edges, log-spaced
        xEdges = np.logspace(
                np.log10(sanX.min()),
                np.log10(sanX.max() + np.diff(sanX)[-1]/100.), #include last point
                nBin + 1)

        # loop over bins:
        for bini in range(nBin):
            fBin[bini], fuBin[bini], x0Bin[bini] = None, None, None # default
            fMask = ((sanX >= xEdges[bini]) & (sanX < xEdges[bini + 1]))
            fInBin, fuInBin = self.f.sanitized[fMask], self.f.sanitizedU[fMask]
            x0InBin = self.x0.sanitized[fMask]
            if fMask.sum() == 1:
                fBin[bini], fuBin[bini], x0Bin[bini] = fInBin, fuInBin, x0InBin
                validMask[bini] = True
            elif fMask.sum() > 1:
                fBin[bini], x0Bin[bini] = fInBin.mean(), x0InBin.mean()
                validMask[bini] = True
                # uncertainties are a bit more elaborate:
                fuBin[bini] = np.maximum(
                        fInBin.std(ddof = 1) / np.sqrt(1. * fMask.sum()), # SEM
                        np.sqrt( (fuInBin**2).sum() / fMask.sum() ) #propagated. unc.
                        )

        # remove empty bins:
        validi = (True - np.isnan(fBin))
        validi[np.argwhere(validMask != True)] = False
        # store values:
        self._fBin = DataVector(u'Ib', fBin[validi], rawU = fuBin[validi], 
                unit = self.f.unit
                ) 
        self._x0Bin = DataVector(u'qb', 
                self.x0.unit.toDisplay(x0Bin[validi]), 
                unit = self.x0.unit
                )
        logging.info("Rebinning procedure completed: {} bins.".format(validi.sum()))

    # These are to be set by the particular application dataset: 
    # i.e.: x = q, y = psi, f = I for SAS, x = tau, f = (G1 - 1) for DLS
    # derived classes may have an alias getter for (x0, f, …)

    @property
    def x0(self):
        """First sampling vector."""
        return self._x0

    @x0.setter
    def x0(self, vec):
        assert vec is None or isinstance(vec, DataVector)
        self._x0 = vec

    @property
    def x1(self):
        """Second sampling vector."""
        return self._x1

    @x1.setter
    def x1(self, vec):
        assert vec is None or isinstance(vec, DataVector)
        self._x1 = vec

    @property
    def x2(self):
        """Third sampling vector."""
        return self._x2

    @x2.setter
    def x2(self, vec):
        assert vec is None or isinstance(vec, DataVector)
        self._x2 = vec

    @property
    def f(self):
        """The measurement vector."""
        return self._f

    @f.setter
    def f(self, vec):
        assert vec is None or isinstance(vec, DataVector)
        self._f = vec
        self._initMask()
        self._propagateMask()

    @property
    def x0Bin(self):
        """Binned variant of First sampling vector."""
        return self._x0Bin

    @x0Bin.setter
    def x0Bin(self, vec):
        logging.error("x0Bin can only be set by the (re-)binning method")

    @property
    def fBin(self):
        """Binned variant of measurement vector."""
        return self._fBin

    @fBin.setter
    def fBin(self, vec):
        logging.error("fBin can only be set by the (re-)binning method")
    # other common meta data

    @classproperty
    @classmethod
    def sourceName(cls):
        """Returns the name of the measurement method."""
        raise NotImplementedError

    @property
    def filename(self):
        return self._filename

    def setFilename(self, fn):
        """Stores the absolute path to this data file.
        Should be reviewed when data sets can be created from several files."""
        if fn is None or not os.path.isfile(fn):
            return
        self._filename = os.path.abspath(fn)

    @abstractproperty
    def count(self):
        raise NotImplementedError

    @property
    def is2d(self):
        """Returns true if this dataset contains two-dimensional data with
        psi information available."""
        res = False
        try:
            res = isinstance(self.x1, DataVector)
        except NotImplementedError:
            pass
        return res

    def accumulate(self, others):
        return None

    @property
    def config(self):
        return self._config

    def setConfig(self, config = None):
        """Set the extended configuration for this data and returns true if
        the configuration was different and an update was necessary."""
        if not isinstance(config, self.configType):
            return False
        # always replacing the config if it's valid, it originates from here
        self._config = config
        self.updateConfig()
        return True

    def updateConfig(self):
        """Updates the config object based on this data set. All callbacks are
        run right after this method in setConfig()."""
        self.config.is2d = self.is2d # forward if we are 2d or not
        self.config.register("x0limits", self._onLimitsUpdate)
        self.config.register("x0Clipping", self._onClippingUpdate)
        self.config.register("fMasks", self._onFMasksUpdate)
        descr = self.config.x0Low.displayName().format(x0 = self.x0.name)
        self.config.x0Low.setDisplayName(descr)
        descr = self.config.x0LowClip.displayName().format(x0 = self.x0.name)
        self.config.x0LowClip.setDisplayName(descr)
        descr = self.config.x0High.displayName().format(x0 = self.x0.name)
        self.config.x0High.setDisplayName(descr)
        descr = self.config.fMaskZero.displayName().format(f = self.f.name)
        self.config.fMaskZero.setDisplayName(descr)
        descr = self.config.fMaskNeg.displayName().format(f = self.f.name)
        self.config.fMaskNeg.setDisplayName(descr)
        # FIXME: Problem with a many2one relation (many data sets, one config)
        #        -> What is the valid range supposed to be?
        #           Atm, the smallest common range wins. [ingo]
        self.config.setX0ValueRange(
                (self.x0.siData.min(), self.x0.siData.max()))
        self._excludeInvalidX0()
        self.reBin()
        if not self.is2d:
            return # self.x1 will be None
        self.config.register("x1limits", self._onLimitsUpdate)
        descr = self.config.x1Low.displayName().format(x1 = self.x1.name)
        self.config.x1Low.setDisplayName(descr)
        descr = self.config.x1High.displayName().format(x1 = self.x1.name)
        self.config.x1High.setDisplayName(descr)
        self.config.setX1ValueRange(
                (self.x1.siData.min(), self.x1.siData.max()))

    def _excludeInvalidX0(self):
        validX0Idx = 0 # get the first data point index above 0
        while self.x0.siData[validX0Idx] <= 0.0:
            validX0Idx += 1
        if self.config.x0LowClip() < validX0Idx:
            self.config.x0LowClip.setValue(validX0Idx)

    @abstractproperty
    def configType(self):
        """Returns a compatible DataConfig type."""
        raise NotImplementedError

    @abstractproperty
    def modelType(self):
        """Returns a compatible ScatteringModel type."""
        raise NotImplementedError

    def _initMask(self):
        # init indices: index array is more flexible than boolean masks
        if self.f is None:
            return
        self._validMask = np.isfinite(self.f.siData)

    def _propagateMask(self):
        # store
        validIndices = np.argwhere(self._validMask)[:,0]
        # pass on all valid indices to the parameters
        self.f.validIndices = validIndices
        self.x0.validIndices = validIndices
        if isinstance(self.x1, DataVector):
            self.x1.validIndices = validIndices
        if isinstance(self.x2, DataVector):
            self.x2.validIndices = validIndices
        # add onMaskUpdate() or validIndicesUpdated() callback here

    def _applyFMasks(self):
        # Optional masking of negative intensity
        if self.config.fMaskZero():
            # FIXME: compare with machine precision (EPS)?
            self._validMask &= (self.f.siData != 0.0)
        if self.config.fMaskNeg():
            self._validMask &= (self.f.siData > 0.0)

    def _applyClipping(self):
        # apply explicit clipping first
        end = max(0, self.config.x0LowClip())
        self._validMask[0:end] = False

    def _applyLimits(self):
        # clip to q bounds
        self._validMask &= (self.x0.siData >= self.config.x0Low())
        self._validMask &= (self.x0.siData <= self.config.x0High())
        # clip to psi bounds
        if not self.is2d:
            return
        # -> is it important to use '>' here, instead of '>=' for x0?
        self._validMask &= (self.x1.siData >  self.config.x1Low())
        self._validMask &= (self.x1.siData <= self.config.x1High())

    def _onFMasksUpdate(self, *args):
        self._initMask()
        self._applyFMasks()
        self._applyLimits()
        self._propagateMask()

    def _onLimitsUpdate(self, *args):
        self._initMask()
        self._applyFMasks()
        self._applyLimits()
        self._propagateMask()
        # update dependent configuration values
        if len(self.x0.validIndices):
            self.config.x0LowClip.setValue(self.x0.validIndices.min())

    def _onClippingUpdate(self, *args):
        self._initMask()
        self._applyFMasks()
        self._applyClipping()
        self._propagateMask()
        # Parameter.setValue() calls the limits callback only once:
        #  It stops calling back in _onLimitsUpdate() because the clipping
        #  value does not change further (no update needed)
        # -> vice versa at the end of _onLimitsUpdate() above
        self.config.x0Low.setValue(self.x0.sanitized.min())

    def __init__(self, **kwargs):
        super(DataObj, self).__init__(**kwargs)

    def __eq__(self, other):
        return (np_all(self.rawArray == other.rawArray)
                and self.title == other.title
                and self.filename == other.filename)

    def __neq__(self, other):
        return not self.__eq__(other)

if __name__ == "__main__":
    import doctest
    doctest.testmod()

# vim: set ts=4 sts=4 sw=4 tw=0:
