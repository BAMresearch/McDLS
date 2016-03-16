# -*- coding: utf-8 -*-
# dataobj/sasconfig.py

from __future__ import absolute_import # PEP328

from abc import ABCMeta, abstractmethod, abstractproperty
import numpy as np
from bases.algorithm import AlgorithmBase
from utils.parameter import Parameter
from utils.units import (ScatteringIntensity, ScatteringVector, Angle,
                         Fraction, NoUnit)
from utils import clip
from dataobj import DataConfig
import logging

class SmearingConfig(AlgorithmBase):
    """Abstract base class, can't be instantiated."""
    __metaclass__ = ABCMeta
    _qOffset = None # integration point positions, depends on beam profile
    _weights = None # integration weight per position, depends on beam profile
    locs = None # integration location matrix, depends on collType
    shortName = "SAS smearing configuration"
    parameters = (
        Parameter("doSmear", False, unit = NoUnit(),
            displayName = "Apply smearing correction",
            ),
        Parameter("nSteps", 25, unit = NoUnit(),
            displayName = "number of smearing points around each q",
            valueRange = (0, 1000)),
        # 2-d collimated systems require a different smearing than slit-collimated data
        Parameter("twoDColl", False, unit = NoUnit(),
            displayName = "Slit-smeared data (unchecked), or 2D-averaged data (checked)",
            ),
#        Parameter("collType", u"Slit", unit = NoUnit(),
#            displayName = "Type of collimation leading to smearing",
#            valueRange = [u"Slit", u"Pinhole", u"Rectangular", u"None"])
    )

    @abstractmethod
    def updateQUnit(self, newUnit):
        pass

    @abstractmethod
    def updateQLimits(self, qLimit):
        pass

    @property
    def qOffset(self):
        return self._qOffset

    @property
    def weights(self):
        return self._weights

    @property
    def prepared(self):
        return self._qOffset, self._weights

    def copy(self):
        other = super(SmearingConfig, self).copy()
        if self.qOffset is not None:
            other._qOffset = self._qOffset.copy()
        if self.weights is not None:
            other._weights = self._weights.copy()
        return other

    def __str__(self):
        s = [str(id(self)) + " " + super(SmearingConfig, self).__str__()]
        s.append("  qOffset: {}".format(self.qOffset))
        s.append("  weights: {}".format(self.weights))
        return "\n".join(s)

# SmearingConfig.factory() # not sure if required

class TrapezoidSmearing(SmearingConfig):
    parameters = (
        Parameter("Umbra", 0., unit = NoUnit(), # unit set outside
            displayName = "top width of <br />trapezoidal beam profile",
            description = "full top width of the trapezoidal beam profile (horizontal for slit-collimated systems, circularly averaged for 2D pinhole and rectangular slit)",
            valueRange = (0., np.inf), decimals = 1),
        Parameter("Penumbra", 0., unit = NoUnit(), # unit set outside
            displayName = "bottom width of <br />trapezoidal beam profile",
            description = "full bottom width of the trapezoidal beam profile horizontal for slit-collimated systems, circularly averaged for 2D pinhole and rectangular slit)",
            valueRange = (0., np.inf), decimals = 1),
    )

    @property
    def showParams(self):
        lst = ["Umbra", "Penumbra"]
        return lst + [name
                for name in super(TrapezoidSmearing, self).showParams
                    if name not in lst]

    def halfTrapzPDF(self, x, c, d):
        # this trapezoidal PDF is only defined from X >= 0, and is assumed
        # to be mirrored around that point. 
        # Note that the integral of this PDF from X>0 will be 0.5. 
        # source: van Dorp and Kotz, Metrika 2003, eq (1) 
        # using a = -d, b = -c
        logging.debug("halfTrapzPDF called")
        assert(d > 0.)
        x = abs(x)
        pdf = x * 0.
        pdf[x < c] = 1.
        if d > c:
            pdf[(c <= x) & (x < d)] = (1./(d - c)) * (d - x[(c <= x) & (x < d)])
        norm = 1./(d + c)
        pdf *= norm
        return pdf, norm

    def setIntPoints(self, q):
        """ sets smearing profile integration points for trapezoidal slit. 
        Top (umbra) of trapezoid has full width xt, bottom of trapezoid 
        (penumbra) has full width.
        Since the smearing function is assumed to be symmetrical, the 
        integration parameters are calculated in the interval [0, xb/2]
        """
        n, xt, xb = self.nSteps(), self.Umbra(), self.Penumbra()
        logging.debug("setIntPoints called with n = {}".format(n))

        # following qOffset is used for Pinhole and Rectangular
        qOffset = np.logspace(np.log10(q.min() / 5.),
                np.log10(xb / 2.), num = np.ceil(n / 2.))
        qOffset = np.concatenate((-qOffset[::-1], [0,], qOffset)) 
        if not self.twoDColl():
            # overwrite prepared integration steps qOffset:
            qOffset = np.logspace(np.log10(q.min() / 5.),
                    np.log10(xb / 2.), num = n)
            # tack on a zero at the beginning
            qOffset = np.concatenate(([0,], qOffset)) 

        y, dummy = self.halfTrapzPDF(qOffset, xt, xb)

        # volume fraction still off by a factor of two (I think). Can be 
        # fixed by multiplying y with 0.5, but need to find it first in eqns. 
                # volume fraction still off by a factor of two (I think). Can be 
                # fixed by multiplying y with 0.5, but need to find it first in eqns. 
        self._qOffset, self._weights = qOffset, y

    def updateQUnit(self, newUnit):
        assert isinstance(newUnit, ScatteringVector)
        self.Umbra.setUnit(newUnit)
        self.Penumbra.setUnit(newUnit)

    def updatePUnit(self, newUnit):
        assert isinstance(newUnit, Angle)
        # TODO

    def updateQLimits(self, qLimit):
        qLow, qHigh = qLimit
        self.Umbra.setValueRange((0., 2. * qHigh))
        self.Penumbra.setValueRange((0., 2. * qHigh))

    def updateSmearingLimits(self, q):
        qHigh = q.max()
        lowLim = diff(q).min()
        self.Umbra.setValueRange((lowLim, 2. * qHigh))
        self.Penumbra.setValueRange((lowLim, 2. * qHigh))

    def updatePLimits(self, pLimit):
        pLow, pHigh = pLimit
        # TODO

    def __init__(self):
        super(TrapezoidSmearing, self).__init__()
        self.Umbra.setOnValueUpdate(self.onUmbraUpdate)

    def onUmbraUpdate(self):
        """Value in umbra will not exceed available q."""
        # value in Penumbra must not be smaller than Umbra
        self.Penumbra.setValueRange((self.Umbra(), self.Penumbra.max()))

TrapezoidSmearing.factory()

class SASConfig(DataConfig):
    # TODO: fix UI elsewhere for unit selection along to each input and forward
    #       that to the DataVector
    _iUnit = NoUnit()
    _qUnit = NoUnit()
    _pUnit = NoUnit()
    _smearing = None
    shortName = "SAS data configuration"

    parameters = (
        Parameter("eMin", Fraction(u"%").toSi(1.), unit = Fraction(u"%"),
            displayName = "minimum uncertainty estimate",
            valueRange = (0., 1.), decimals = 1),
    )

    @property
    def callbackSlots(self):
        return super(SASConfig, self).callbackSlots | set((
            "qunit", "punit", "iunit", "eMin"))

    def updateEMin(self):
        self.callback("eMin", self.eMin())

    def setX0ValueRange(self, limit):
        """Sets available range of loaded data."""
        super(SASConfig, self).setX0ValueRange(limit)
        if self.smearing is None:
            return
        self.smearing.updateQLimits((self.x0Low(), self.x0High())) # correct?
        # self.smearing.updateSmearingLimits(self.q) # correct?

    def setX1ValueRange(self, limit):
        super(SASConfig, self).setX1ValueRange(limit)
        # TODO

    @property
    def iUnit(self):
        return self._iUnit

    @property
    def qUnit(self):
        return self._qUnit

    @property
    def pUnit(self):
        return self._pUnit

    @iUnit.setter
    def iUnit(self, newUnit):
        assert isinstance(newUnit, ScatteringIntensity)
        self._iUnit = newUnit
        self.callback("iunit", newUnit)

    @qUnit.setter
    def qUnit(self, newUnit):
        assert isinstance(newUnit, ScatteringVector)
        self._qUnit = newUnit
        self.callback("qunit", newUnit)

    @pUnit.setter
    def pUnit(self, newUnit):
        assert isinstance(newUnit, Angle)
        self._pUnit = newUnit
        self.callback("punit", newUnit)

    @property
    def smearing(self):
        return self._smearing

    @smearing.setter
    def smearing(self, newSmearing):
        assert isinstance(newSmearing, SmearingConfig)
        self._smearing = newSmearing

    def prepareSmearing(self, q):

        assert( isinstance(q, np.ndarray))
        assert( q.ndim == 1)

        if self.smearing is None:
            return q
        if self.smearing.Penumbra() == 0.:
            return q
        if not self.smearing.doSmear():
            logging.debug("smearing disabled")
            return q
        self.smearing.setIntPoints(q)
        qOffset, weights = self.smearing.prepared
        #print >>sys.__stderr__, "prepareSmearing"
        #print >>sys.__stderr__, unicode(self)
        # calculate the intensities at sqrt(q**2 + qOffset **2)
        if not self.smearing.twoDColl(): # slit collimation
            logging.debug("prepareSmearing called for slit collimation")
            logging.debug("q.shape: {}, qOffset.shape: {}".format(q.shape, qOffset.shape))
            return np.sqrt(np.add.outer(q **2, qOffset **2))
        else: 
            logging.debug("prepareSmearing called for pinhole collimation")
            # Non-slit-smeared instruments, using azimuthally averaged
            # 2D-pattern (assumed!) with equally averaged beam profile.
            logging.debug("q.shape: {}, qOffset.shape: {}".format(q.shape, qOffset.shape))
            logging.debug("qOffset.min: {}, qOffset.max: {}"
                    .format(qOffset.min(), qOffset.max()))
            return np.add.outer(q, qOffset)

    def copy(self):
        other = super(SASConfig, self).copy(smearing = self.smearing.copy())
        other.iUnit = self.iUnit
        other.qUnit = self.qUnit
        other.pUnit = self.pUnit
        return other

    def __init__(self, *args, **kwargs):
        super(SASConfig, self).__init__()
        smearing = kwargs.pop("smearing", None)
        if smearing is None:
            smearing = TrapezoidSmearing()
        self.smearing = smearing
        self.register("qunit", self.x0Low.setUnit)
        self.register("qunit", self.x0High.setUnit)
        self.register("punit", self.x1Low.setUnit)
        self.register("punit", self.x1High.setUnit)
        if self.smearing is not None:
            self.register("qunit", self.smearing.updateQUnit)
            self.register("punit", self.smearing.updatePUnit)
            self.register("x0limits", self.smearing.updateQLimits)
            # self.register("x0limits", self.smearing.updateSmearingLimits)
            self.register("x1limits", self.smearing.updatePLimits)
        self.iUnit = ScatteringIntensity(u"(m sr)⁻¹")
        self.qUnit = ScatteringVector(u"nm⁻¹")
        self.pUnit = Angle(u"°")
        self.eMin.setOnValueUpdate(self.updateEMin)

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        equal = ((self.smearing == other.smearing) and
                super(SASConfig, self).__eq__(other))
        return equal

    def __str__(self):
        return "\n".join((
            super(SASConfig, self).__str__(),
            unicode(self.smearing)
        ))

SASConfig.factory() # check class variables

# vim: set ts=4 sts=4 sw=4 tw=0:
