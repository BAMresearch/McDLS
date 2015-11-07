# -*- coding: utf-8 -*-
# dataobj/sasdata.py

"""
Represents data associated with a measurement by small angle scattering (SAS).
Some examples and tests.

>>> import numpy
>>> testdata = numpy.random.rand(4,4)
>>> testtitle = "some title"
>>> from sasdata import SASData

Testing
>>> first = SASData(testtitle, testdata)
>>> first.title == testtitle
True
>>> numpy.all(first.rawArray == testdata)
True
"""

from __future__ import absolute_import # PEP328
import logging
import numpy as np # For arrays
from utils import classproperty
from utils.units import Length, ScatteringVector, ScatteringIntensity, Angle
from dataobj import DataObj, SASConfig

class SASData(DataObj):
    """Represents one set of data from a unique source (a file, for example).
    """
    _sizeEst = None
    _shannonChannelEst = None
    _eMin = 0.01 # minimum possible error (1%)
    _uncertainty = None
    _rUnit = None # defines units for r used in sizeest
    _qClipRange = [-np.inf, np.inf] # Q clip range for this dataset
    _pClipRange = [-np.inf, np.inf] # psi clip range for this dataset
    _maskZeroInt = False # mask zero intensity values (updated by mcsasparam)
    _maskNegativeInt = False # mask negative intensity values (updated by mcsasparam)
    _validIndices = None # valid indices generated by the prepareClipMask function
    _config = None

    # for smearing:
    _doSmear = False # smear the data
    _slitUmbra = None # top width of the trapezoidal beam length profile
    _slitPenumbra = None # top width of the trapezoidal beam length profile
    slitShape = "trapezoid" # default
    locs = None # substitute array of Q-values to calculate for smearing
    dU = None # array of smear integration values
    weightFunc = None # array of smearing weighting values.
    nSmearSteps = 25 # number of smearing points

    # intensity related helpers

    @property
    def i(self):
        """Measured intensity at q."""
        return self.iOrigin.copy()[self.validIndices]

    @property
    def iOrigin(self):
        return self.iUnit.toSi(self.rawArray[:, 1])

    @property
    def iUnit(self):
        return self.config.iUnit

    # scattering vector

    @property
    def q(self):
        """Q-Vector at which the intensities are measured."""
        return self.qOrigin.copy()[self.validIndices]

    @property
    def qOrigin(self):
        return self.qUnit.toSi(self.rawArray[:, 0])

    @property
    def qUnit(self):
        return self.config.qUnit

    @property
    def qMin(self):
        """Returns minimum q from data or qClipRange, whichever is larger."""
        return np.maximum(self.qClipRange[0], 
                self.qOrigin.min())

    @qMin.setter
    def qMin(self, newParam):
        """Value in cliprange will not exceed available q."""
        self._qClipRange[0] = np.maximum(
                newParam, self.qOrigin.min())
        self._prepareValidIndices()

    @property
    def qMax(self):
        return np.minimum(self.qClipRange[1], 
                self.qOrigin.max())

    @qMax.setter
    def qMax(self, newParam):
        """Value in cliprange will not exceed available q."""
        self._qClipRange[1] = np.minimum(
               newParam, self.qOrigin.max())
        self._prepareValidIndices()

    @property
    def qClipRange(self):
        return self._qClipRange

    @qClipRange.setter
    def qClipRange(self, newParam):
        """Expects a range in self.qUnits, not SI units!"""
        if not (np.size(newParam) == 2):
            logging.error('qClipRange must be supplied with two-element vector')
        else:
            self.qMin = np.min(newParam)
            self.qMax = np.max(newParam)

    @property
    def qLimsString(self):
        return u"{0:.3g} ≤ Q ({qMagnitudeName}) ≤ {1:.3g}".format(
                self.qUnit.toDisplay(self.qMin),
                self.qUnit.toDisplay(self.qMax),
                qMagnitudeName = self.qUnit.displayMagnitudeName)

    # uncertainty on the intensities from data file

    @property
    def e(self):
        """Uncertainty or Error of the intensity at q loaded from file."""
        return self.eOrigin.copy()[self.validIndices]

    @property
    def eOrigin(self):
        return self.iUnit.toSi(self.rawArray[:, 2])

    @property
    def eMin(self):
        return self._eMin

    @eMin.setter
    def eMin(self, value):
        value = float(value) 
        assert((value > 0.) and (value < 1.))
        self._eMin = value
        #update uncertainty in case eMin was adjusted
        self._prepareUncertainty() 
        self._prepareValidIndices()

    # sanitized uncertainty on the intensities

    @property
    def u(self):
        """Corrected uncertainty or error of the intensity at q."""
        return self.uOrigin.copy()[self.validIndices]

    @property
    def uOrigin(self):
        return self._uncertainty

    # psi scattering vector for 2D data

    @property
    def p(self):
        """Psi-Vector."""
        return self.rawArray[:, 3]

    @property
    def pUnit(self):
        return self.config.pUnit

    @property
    def pMin(self):
        """Returns minimum psi from data or psiClipRange, whichever is larger."""
        if self.is2d:
            return np.maximum(self.pClipRange[0], 
                    self.pOrigin.min())
        else:
            return self.pClipRange[0]

    @pMin.setter
    def pMin(self, newParam):
        """Value in cliprange will not exceed available psi."""
        if self.is2d:
            self._pClipRange[0] = np.maximum(
                newParam, self.pOrigin.min())
        else:
            self._pClipRange[0] = self.pUnit.toSi(newParam)
        self._prepareValidIndices()

    @property
    def pMax(self):
        if self.is2d:
            return np.minimum(self.pClipRange[1], 
                    self.pOrigin.max())
        else:
            return self.pClipRange[1]

    @pMax.setter
    def pMax(self, newParam):
        """Value in cliprange will not exceed available psi."""
        if self.is2d:
            self._pClipRange[1] = np.minimum(
                newParam, self.pOrigin.max())
        else:
            self._pClipRange[1] = self.pUnit.toSi(newParam)
        self._prepareValidIndices()

    @property
    def pClipRange(self):
        return self._pClipRange

    @pClipRange.setter
    def pClipRange(self, newParam):
        if not (np.size(newParam) == 2):
            logging.error('pClipRange must be supplied with two-element vector')
        else:
            self.pMin(np.min(newParam))
            self.pMax(np.min(newParam))

    @property
    def pLimsString(self):
        return u"{0:.3g} ≤ psi ({psiMagnitudeName}) ≤ {1:.3g}".format(
                self.pUnit.toDisplay(self.pMin),
                self.pUnit.toDisplay(self.pMax),
                pMagnitudeName = self.pUnit.displayMagnitudeName)

    # for smearing, define beam profile (one direction only for now)
    @property
    def doSmear(self):
        return self._doSmear
    
    @doSmear.setter
    def doSmear(self, value):
        self._doSmear = value
        if value and ((self.slitPenumbra > 0.) or (self.slitUmbra > 0.)):
            self._prepSmear()

    @property
    def slitUmbra(self):
        return self._slitUmbra

    @slitUmbra.setter
    def slitUmbra(self, newParam):
        """Value in slitUmbra will not exceed available q."""
        self._slitUmbra = np.minimum(
               newParam, self.qOrigin.max())
        # value in Penumbra must not be smaller than Umbra
        self._slitPenumbra = np.maximum(self._slitUmbra, self._slitPenumbra)
        if self.doSmear:
            self._prepSmear()

    @property
    def slitPenumbra(self):
        return self._slitPenumbra

    @slitPenumbra.setter
    def slitPenumbra(self, newParam):
        """Value in slitPenumbra will not exceed available q, and not be smaller 
        than Umbra."""
        self._slitPenumbra = np.minimum(
               newParam, self.qOrigin.max())
        self._slitPenumbra = np.maximum(
               newParam, self._slitUmbra)
        if self.doSmear:
            self._prepSmear()

    # general information on this data set

    @property
    def count(self):
        return len(self.q)

    @property
    def validIndices(self):
        if self._validIndices is not None:
            return self._validIndices
        else:
            # valid indices not set yet
            self._prepareValidIndices()
            return self._validIndices

    @property
    def is2d(self):
        """Returns true if this dataset contains two-dimensional data with
        psi information available."""
        return self.rawArray.shape[1] > 3 # psi column is present

    @property
    def hasError(self):
        """Returns True if this data set has an error bar for its
        intensities."""
        return self.rawArray.shape[1] > 2

    # general info texts for the UI

    @property
    def dataContent(self):
        """shows the content of the loaded data: Q, I, IErr, etc"""
        content = []
        if self.q is not None:
            content.append('Q')
        if self.i is not None:
            content.append('I')
        if self.hasError:
            content.append('IErr')
        if self.is2d:
            content.append('Psi')
        return ", ".join(content)

    @classproperty
    @classmethod
    def displayDataDescr(cls):
        return ("Filename ", "Data points ", "Data content ", 
                "Q limits ", "Est. sphere size ", "Recommended number of bins ")

    @classproperty
    @classmethod
    def displayData(cls):
        return ("title", "count", "dataContent", 
                "qLimsString", "sphericalSizeEstText", "shannonChannelEstText")

    def sphericalSizeEst(self):
        return self._sizeEst

    @property
    def rUnit(self):
        return self._rUnit

    @property
    def sphericalSizeEstText(self):
        return u"{0:.3g} ≤ R ({rUnitName}) ≤ {1:.3g} ".format(
                *self.rUnit.toDisplay(self.sphericalSizeEst()),
                rUnitName = self.rUnit.displayMagnitudeName)

    def shannonChannelEst(self):
        return int(self._shannonChannelEst)

    @property
    def shannonChannelEstText(self):
        return u"≤ {0:2d} bins ".format(self.shannonChannelEst())

    # masking helpers for sanitizing input data

    @property
    def maskZeroInt(self):
        return bool(self._maskZeroInt)

    @maskZeroInt.setter
    def maskZeroInt(self, value):
        self._maskZeroInt = bool(value)
        self._prepareValidIndices()

    @property
    def maskNegativeInt(self):
        return bool(self._maskNegativeInt)

    @maskNegativeInt.setter
    def maskNegativeInt(self, value):
        self._maskNegativeInt = bool(value)
        self._prepareValidIndices()

    @property
    def config(self):
        return self._config

    def __init__(self, **kwargs):
        super(SASData, self).__init__(**kwargs)
        #set unit definitions for display and internal units
        self._config = SASConfig()
        import sys
        print >>sys.__stderr__, self.__class__.__name__
        print >>sys.__stderr__, unicode(self.config)
        self._rUnit = Length(u"nm")

        self._prepareValidIndices()
        self._prepareSizeEst()
        self._shannonChannelEst = self.q.max() / self.q.min()
        self._prepareUncertainty()

    def _prepareSizeEst(self):
        self._sizeEst = np.pi / np.array([self.q.max(),
                                          abs(self.q.min()) ])

    def _prepareUncertainty(self):
        self._uncertainty = self.eMin * self.iOrigin
        minUncertaintyPercent = self.eMin * 100.
        if not self.hasError:
            logging.warning("No error column provided! Using {}% of intensity."
                            .format(minUncertaintyPercent))
        else:
            count = sum(self._uncertainty > self.eOrigin)
            if count > 0:
                logging.warning("Minimum uncertainty ({}% of intensity) set "
                                "for {} datapoints.".format(
                                minUncertaintyPercent, count))
            self._uncertainty = np.maximum(self._uncertainty, self.eOrigin)
        # reset invalid uncertainties to np.inf
        invInd = (True - np.isfinite(self._uncertainty))
        self._uncertainty[invInd] = np.inf

    def _prepSmear(self):


        # def prepSmear(self, data, slitShape = "trapezoid", shapeParam = [0., 0.], nSmearSteps = 25):

        def squareSlit(q, shapeParam, n):
            """ Testing purposes only, since trapezoid encompasses square slit """

            slitWidth = shapeParam[0]
            # prepare integration steps dU:
            dU = np.logspace(np.log10(q.min() / 10.),
                    np.log10(slitWidth / 2.), num = n - 1)
            dU = np.concatenate(([0,], dU)) [np.newaxis, :]
            return dU, np.ones((dU.size,)) / slitWidth

        def trapzSlit(q, shapeParam, n):
            """ defines integration over trapezoidal slit. Top of trapezoid 
            has width xt, bottom of trapezoid has width xb. Note that xb > xt"""
            xt, xb = shapeParam[0], shapeParam[1]

            # ensure things are what they are supposed to be
            assert (xt >= 0.)
            if xb < xt:
                xb = xt # should use square profile in this case.

            # prepare integration steps dU:
            dU = np.logspace(np.log10(q.min() / 10.),
                    np.log10(xb / 2.), num = n)
            dU = np.concatenate(([0,], dU)) [np.newaxis, :]

            if xb == xt: 
                y = 1. - (dU * 0.)
            else:
                y = 1. - (dU - xt) / (xb - xt)

            y = np.clip(y, 0., 1.)
            y[dU < xt] = 1.
            Area = (xt + 0.5 * (xb - xt))
            return dU, y / Area
        
        # now we do the actual smearing preparation
        assert isinstance(self.q, np.ndarray)
        assert (self.q.ndim == 1)

        # define smearing profile
        if self.slitShape == "trapezoid":
            slitFcn = trapzSlit
        else:
            slitFcn = squareSlit

        self.dU, self.weightFunc = slitFcn(self.q, 
                [self.slitUmbra, self.slitPenumbra], 
                self.nSmearSteps)

        # calculate the intensities at sqrt(q**2 + dU **2)
        self.locs = np.sqrt(np.add.outer(self.q **2, self.dU[0,:]**2)) 
    
    def _prepareValidIndices(self):
        """
        If q and/or psi limits are supplied in the dataset,
        prepares a clipmask for the full dataset
        """
        # init indices: index array is more flexible than boolean masks
        bArr = np.isfinite(self.qOrigin)

        # Optional masking of negative intensity
        if self.maskZeroInt:
            # FIXME: compare with machine precision (EPS)?
            bArr &= (self.iOrigin != 0.0)
        if self.maskNegativeInt:
            bArr &= (self.iOrigin > 0.0)

        # clip to q bounds
        bArr &= (self.qOrigin >= self.qMin)
        bArr &= (self.qOrigin <= self.qMax)
        # clip to psi bounds
        if self.is2d:
            bArr &= (self.pOrigin > self.pMin)
            bArr &= (self.pOrigin <= self.pMax)

        # store
        self._validIndices = np.argwhere(bArr)[:,0]
        self._prepareSizeEst() # recalculate based on limits. 

if __name__ == "__main__":
    import doctest
    doctest.testmod()

# vim: set ts=4 sts=4 sw=4 tw=0:
