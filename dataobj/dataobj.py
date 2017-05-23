# -*- coding: utf-8 -*-
# scattering/dataobj.py

"""
Represents input data associated with a measurement.
"""

from __future__ import absolute_import # PEP328
from builtins import range
import os # Miscellaneous operating system interfaces
from numpy import all as np_all
import numpy as np

from abc import ABCMeta, abstractproperty
from bases.dataset import DataSet, DisplayMixin
from dataobj.datavector import DataVector
from utils import classproperty, hashNumpyArray
import logging
from future.utils import with_metaclass

class DataObj(with_metaclass(ABCMeta, type('NewBase', (DataSet, DisplayMixin), {}))):
    """General container for data loaded from file. It offers specialised
    methods to derive information from the provided data.
    """
    _filename = None
    _config = None
    _validMask = None
    _x0 = None
    _x1 = None
    _x2 = None
    _f  = None

    # The following are to be set by the particular application dataset: 
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

    @classproperty
    @classmethod
    def sourceName(cls):
        """Returns the name of the measurement method."""
        raise NotImplementedError

    @property
    def sampleName(self):
        return ""

    @property
    def filename(self):
        return self._filename

    def setFilename(self, fn):
        """Stores the absolute path to this data file.
        Should be reviewed when data sets can be created from several files."""
        if fn is None or not os.path.isfile(fn):
            return
        self._filename = os.path.abspath(fn)

    @property
    def seriesKey(self):
        """The Name of the DataObj property to use as series key,
        hard-coded for now, assuming it exists.
        It allows to let the user chose from a generated list of properties
        (todo)."""
        return "title"

    @property
    def seriesKeyName(self):
        """Returns the docstring of the property defined by self.seriesKeyProp.
        """
        try:
            return getattr(type(self), self.seriesKey, None).__doc__
        except AttributeError:
            raise
        return ""

    @property
    def seriesKeyValue(self):
        """Returns the value of the property defined by self.seriesKeyProp."""
        return getattr(self, self.seriesKey, None)

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

    @abstractproperty
    def configType(self):
        """Returns a compatible DataConfig type."""
        raise NotImplementedError

    @property
    def config(self):
        return self._config

    def initConfig(self):
        """Initializes a new data configuration and sets the sample name which
        is used to differentiate different data objects of the same type later
        on."""
        config = self.configType()
        # important to pass the check in setConfig()
        config.sampleName = self.sampleName
        self.setConfig(config)

    def setConfig(self, config = None):
        """Set the configuration of this data object if the type matches."""
        if not isinstance(config, self.configType):
            return # ignore configurations of other types
        if self.config is not None and self.config.sampleName != config.sampleName:
            return # ignore data configurations of other samples
        # always replacing the config if it's valid, it originates from here
        self._config = config
        self.updateConfig()

    def updateConfig(self):
        """Updates the config object based on this data set. All callbacks are
        run right after this method in setConfig()."""
#        self.config.sampleName = self.sampleName # moved to initConfig()
        self.config.is2d = self.is2d # forward if we are 2d or not
        self.config.register("x0limits", self._onLimitsUpdate)
        self.config.register("fMasks", self._onFMasksUpdate)
        self.config.register("fuMin", self._prepareUncertainty)
        self.config.x0Low.formatDisplayName(x0 = self.x0.name)
        self.config.x0High.formatDisplayName(x0 = self.x0.name)
        self.config.updateX0Unit(self.x0.unit)
        self.config.fMaskZero.formatDisplayName(f = self.f.name)
        self.config.fMaskNeg.formatDisplayName(f = self.f.name)
        self.config.fuRel.formatDisplayName(f = self.f.name)
        # FIXME: Problem with a many2one relation (many data sets, one config)
        #        -> What is the valid range supposed to be?
        #           Atm, the smallest common range wins. [ingo]
        self.config.onUpdatedX0(self.x0.siData)
        self._excludeInvalidX0()
        self._reBin()
        if not self.is2d:
            return # self.x1 will be None
        self.config.register("x1limits", self._onLimitsUpdate)
        self.config.x1Low.formatDisplayName(x1 = self.x1.name)
        self.config.x1High.formatDisplayName(x1 = self.x1.name)
        self.config.updateX1Unit(self.x1.unit)
        self.config.onUpdatedX1(self.x1.siData)

    def hdfWrite(self, hdf):
        hdf.writeMember(self, "filename")
        hdf.writeMembers(self, "f", "x0", "x1", "config")

    def _excludeInvalidX0(self):
        validX0Idx = 0 # get the first data point index above 0
        while self.x0.siData[validX0Idx] <= 0.0:
            validX0Idx += 1

    def _prepareUncertainty(self, *dummy):
        """Modifies the uncertainty of the whole range of measured data to be
        above a previously set minimum threshold *fuMin*."""
        siDataU = self.f.unit.toSi(self.f.rawDataU)
        if self.config.fuOne(): # set uncertainties == 1
            siDataU = np.ones_like(self.f.siDataU) * 0.1
            logging.info("Setting uncertainties = {0:.2g}.".format(siDataU.mean()))
        if self.config.fuRel(): # divide the uncertainties by the signal
            logging.info("Dividing uncertainties by the measured signal.")
            siDataU = siDataU / self.f.siData
        # make sure, uncertainties meet the specified minimum
        minUncertaintyPercent = self.config.fuMin() * 100.
        siDataUMin = self.config.fuMin() * self.f.siData
        if not self.hasUncertainties:
            logging.warning("No error column provided! Using {}% of intensity."
                            .format(minUncertaintyPercent))
            self.f.siDataU = siDataUMin
        else:
            upd = np.maximum(siDataU, siDataUMin)
            count = sum(upd <= siDataUMin)
            if count > 0:
                logging.warning("Minimum uncertainty of {}% intensity set "
                                "for {} data points.".format(
                                minUncertaintyPercent, count))
            else:
                logging.info("No data point falls behind minimum uncertainty "
                         "of {}% intensity.".format(minUncertaintyPercent))
            self.f.siDataU = upd
        # reset invalid uncertainties to np.inf
        invInd = (True - np.isfinite(self.f.siDataU))
        self.f.siDataU[invInd] = np.inf

    @property
    def hasUncertainties(self):
        """Returns True if this data set has an error bar for its
        intensities."""
        return self.f.rawDataU is not None or all(self.f.rawDataU == 0.)

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

    def _reBin(self):
        """Rebinning method, to be run (f.ex.) upon every "Start" buttonpress.
        For now, this will rebin using the x0 vector as a base, although the
        binning vector can theoretically be chosen freely.
        """
        sanX = self.x0.sanitized
        if not len(sanX):
            return

        nBin = self.config.nBin.value()
        if not(nBin > 0):
            self.x0.binnedData = None # reset to none if set
            self.f.binnedData = None
            self.f.binnedDataU = None
            return # no need to do the actual rebinning. values stay None.

        logging.info("Initiating binning procedure for {} bins".format(nBin))

        # self._binned = DataVector() once binning finishes.. dataVector can be set once.
        x0Bin = np.zeros(nBin)
        fBin  = np.zeros(nBin)
        fuBin  = np.zeros(nBin)
        validMask = np.zeros(nBin, dtype = bool) #default false

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
            fInBinDDoF = 0
            if len(fInBin) > 1: # prevent division by zero in numpy.std()
                fInBinDDoF = 1
            x0InBin = self.x0.sanitized[fMask]
            if fMask.sum() == 1:
                fBin[bini], fuBin[bini], x0Bin[bini] = fInBin, fuInBin, x0InBin
                validMask[bini] = True
            elif fMask.sum() > 1:
                fBin[bini], x0Bin[bini] = fInBin.mean(), x0InBin.mean()
                validMask[bini] = True
                # uncertainties are a bit more elaborate:
                fuBin[bini] = np.maximum(
                        fInBin.std(ddof = fInBinDDoF) / np.sqrt(1. * fMask.sum()), # SEM
                        np.sqrt( (fuInBin**2).sum() / fMask.sum() ) #propagated. unc.
                        )

        # remove empty bins:
        validi = (True - np.isnan(fBin))
        validi[np.argwhere(validMask != True)] = False
        # store values:
        self.f.binnedData, self.f.binnedDataU = fBin[validi], fuBin[validi]
        self.x0.binnedData = x0Bin[validi] # self.x0.unit.toDisplay(x0Bin[validi])
        logging.info("Rebinning procedure completed: {} bins.".format(validi.sum()))


    def __init__(self, **kwargs):
        super(DataObj, self).__init__(**kwargs)

    def __hash__(self):
        value = hash(self.title) ^ hash(self.filename)
        value ^= hashNumpyArray(self.rawArray)
        return value

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __neq__(self, other):
        return not self.__eq__(other)

if __name__ == "__main__":
    import doctest
    doctest.testmod()

# vim: set ts=4 sts=4 sw=4 tw=0:
