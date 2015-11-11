# -*- coding: utf-8 -*-
# dataobj/sasconfig.py

from __future__ import absolute_import # PEP328

from abc import ABCMeta, abstractmethod, abstractproperty
import numpy
from bases.algorithm import AlgorithmBase
from utils.parameter import Parameter
from utils.mixedmethod import mixedmethod
from utils.units import ScatteringIntensity, ScatteringVector, Angle, NoUnit
from utils import clip

class SmearingConfig(AlgorithmBase):
    """Abstract base class, can't be instantiated."""
    __metaclass__ = ABCMeta
    _dU = None # integration point positions, depends on beam profile
    _weights = None # integration weight per position, depends on beam profile
    locs = None # integration location matrix, depends on collType
    parameters = (
        # not sure if this is the right place: is the nsteps parameter useful
        # for all possible smearing settings? BRP: yes, I think so... 
        Parameter("nSteps", 25, unit = NoUnit(),
            displayName = "number of smearing points around each q",
            valueRange = (0, 1000)),
#        Parameter("collType", u"Slit", unit = NoUnit(),
#            displayName = "Type of collimation leading to smearing",
#            valueRange = [u"Slit", u"Pinhole", u"Rectangular", u"None"])
    )

    @abstractmethod
    def updateQUnit(self, newUnit):
        pass

    @abstractmethod
    def updateQLimits(self, qLow, qHigh):
        pass

    @abstractmethod
    def integrate(self, q):
        # now we do the actual smearing preparation
        assert isinstance(q, numpy.ndarray)
        assert (q.ndim == 1)

    @property
    def dU(self):
        return self._dU

    @property
    def weights(self):
        return self._weights

    @property
    def prepared(self):
        return self._dU, self._weights

    def copy(self):
        other = super(SmearingConfig, self).copy()
        if self.dU is not None:
            other._dU = self._dU.copy()
        if self.weights is not None:
            other._weights = self._weights.copy()
        return other

    def __str__(self):
        s = [str(id(self)) + " " + super(SmearingConfig, self).__str__()]
        s.append("  dU: {}".format(self.dU))
        s.append("  weights: {}".format(self.weights))
        return "\n".join(s)

# SmearingConfig.factory() # not sure if required

import sys

class TrapezoidSmearing(SmearingConfig):
    parameters = (
        Parameter("umbra", 1e9, unit = NoUnit(), # unit set outside
            displayName = "top width of the trapezoidal beam length profile",
            valueRange = (0., numpy.inf), decimals = 1),
        Parameter("penumbra", 2e9, unit = NoUnit(), # unit set outside
            displayName = "bottom width of the trapezoidal beam length profile",
            valueRange = (0., numpy.inf), decimals = 1),
    )

    def _prepSmear(self, q):
        """ prepares the smearing profile for a given collimation configuration. 
        This is supposed to be in the SmearingConfig class """

        # make sure we're getting a valid dataset:
        assert( isinstance(q, np.ndarray))
        assert( q.ndim == 1)

        # prepare the smearing profile
        self.setIntPoints(q)

        if self.collType == u"Slit":
            self.locs = np.sqrt(np.add.outer(q **2, self.dU[0,:] **2))
        elif self.collType == u"None":
            pass
        else:
            raise NotImplementedError

    def setIntPoints(self, q):
        """ sets smearing profile integration points for trapezoidal slit. 
        Top (umbra) of trapezoid has full width xt, bottom of trapezoid 
        (penumbra) has full width.
        Since the smearing function is assumed to be symmetrical, the 
        integration parameters are calculated in the interval [0, xb/2]
        """
        xt, xb = self.umbra, self.penumbra

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
        self._dU, self._weights = dU, (y / Area)
    
    def updateQUnit(self, newUnit):
        assert isinstance(newUnit, ScatteringVector)
        self.umbra.setUnit(newUnit)
        self.penumbra.setUnit(newUnit)

    def updateQLimits(self, qLow, qHigh):
        self.umbra.setValueRange((qLow, qHigh))
        self.penumbra.setValueRange((self.umbra(), qHigh))

    def __init__(self):
        super(TrapezoidSmearing, self).__init__()
        self.umbra.setOnValueUpdate(self.onUmbraUpdate)

    def onUmbraUpdate(self):
        """Value in umbra will not exceed available q."""
        # value in Penumbra must not be smaller than Umbra
        self.penumbra.setValueRange((self.umbra(), self.penumbra.max()))

    def integrate(self, q):
        """ defines integration over trapezoidal slit. Top of trapezoid 
        has width xt, bottom of trapezoid has width xb. Note that xb > xt"""
        super(TrapezoidSmearing, self).integrate(q)
        n, xt, xb = self.nSteps(), self.umbra(), self.penumbra()
        print >>sys.__stderr__, "integrate", xb, xt

        # ensure things are what they are supposed to be
        assert (xt >= 0.)
        if xb < xt:
            xb = xt # should use square profile in this case.

        # prepare integration steps dU:
        dU = numpy.logspace(numpy.log10(q.min() / 10.),
                            numpy.log10(xb / 2.), num = n)
        dU = numpy.concatenate(([0,], dU)) [numpy.newaxis, :]

        if xb == xt: 
            y = 1. - (dU * 0.)
        else:
            y = 1. - (dU - xt) / (xb - xt)

        y = numpy.clip(y, 0., 1.)
        y[dU < xt] = 1.
        area = (xt + 0.5 * (xb - xt))
        self._dU, self._weights = dU, y / area

TrapezoidSmearing.factory()

class SASConfig(AlgorithmBase):
    # set units already for use in parameter definitions and UI
    _iUnit = NoUnit()
    _qUnit = NoUnit()
    _pUnit = NoUnit()
    _smearing = None
    parameters = (
        Parameter("qLow", 0., unit = NoUnit(),
            displayName = "lower q cut-off",
            valueRange = (0., numpy.inf), decimals = 1),
        Parameter("qHigh", numpy.inf, unit = NoUnit(),
            displayName = "upper q cut-off",
            valueRange = (0., numpy.inf), decimals = 1),
    )

    @mixedmethod
    def updateConstraints(self, *args):
        if not self.qLow() < self.qHigh():
            temp = self.qLow()
            self.qLow.setValue(self.qHigh())
            self.qHigh.setValue(temp)
        if self.smearing is None:
            return
        self.smearing.updateQLimits(self.qLow(), self.qHigh())

    # not used yet, need a signal to get the available q-range of data from filelist to datawidget
    def setQRange(self, qMin, qMax):
        """Sets available range of loaded data."""
        self.qLow.setValueRange((qMin, qMax))
        self.qHigh.setValueRange((qMin, qMax))
        if self.smearing is None:
            return
        self.smearing.updateQLimits(self.qLow(), self.qHigh())

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

    @qUnit.setter
    def qUnit(self, newUnit):
        assert isinstance(newUnit, ScatteringVector)
        self._qUnit = newUnit
        self.qLow.setUnit(newUnit)
        self.qHigh.setUnit(newUnit)
        if self.smearing is None:
            return
        self.smearing.updateQUnit(newUnit)

    @pUnit.setter
    def pUnit(self, newUnit):
        assert isinstance(newUnit, Angle)
        self._pUnit = newUnit

    @property
    def smearing(self):
        return self._smearing

    @smearing.setter
    def smearing(self, newSmearing):
        assert isinstance(newSmearing, SmearingConfig)
        self._smearing = newSmearing

    def prepareSmearing(self, q):
        if self.smearing is None:
            return
        self.smearing.integrate(q)
        dU, weights = self.smearing.prepared
        print >>sys.__stderr__, "prepareSmearing"
        print >>sys.__stderr__, unicode(self)
        # calculate the intensities at sqrt(q**2 + dU **2)
        return numpy.sqrt(numpy.add.outer(q**2, dU[0,:]**2))

    def copy(self):
        other = super(SASConfig, self).copy()
        other.iUnit = self.iUnit
        other.qUnit = self.qUnit
        other.pUnit = self.pUnit
        other.smearing = self.smearing.copy()
        return other

    def __init__(self):
        super(SASConfig, self).__init__()
        self.smearing = TrapezoidSmearing()
        self.iUnit = ScatteringIntensity(u"(m sr)⁻¹")
        self.qUnit = ScatteringVector(u"nm⁻¹")
        self.pUnit = Angle(u"°")
        self.qLow.setOnValueUpdate(self.updateConstraints)
        self.qHigh.setOnValueUpdate(self.updateConstraints)

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
