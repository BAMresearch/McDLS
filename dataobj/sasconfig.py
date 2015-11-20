# -*- coding: utf-8 -*-
# dataobj/sasconfig.py

from __future__ import absolute_import # PEP328

from abc import ABCMeta, abstractmethod, abstractproperty
import numpy
from bases.algorithm import AlgorithmBase
from utils.parameter import Parameter
from utils.mixedmethod import mixedmethod
from utils.units import ScatteringIntensity, ScatteringVector, Angle, NoUnit
from utils import clip, isCallable

class SmearingConfig(AlgorithmBase):
    """Abstract base class, can't be instantiated."""
    __metaclass__ = ABCMeta
    _qOffset = None # integration point positions, depends on beam profile
    _weights = None # integration weight per position, depends on beam profile
    locs = None # integration location matrix, depends on collType
    shortName = "SAS smearing configuration"
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
    def updateQLimits(self, qLimit):
        pass

    @abstractmethod
    def integrate(self, q):
        # now we do the actual smearing preparation
        assert isinstance(q, numpy.ndarray)
        assert (q.ndim == 1)

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
            self.locs = np.sqrt(np.add.outer(q **2, self.qOffset[0,:] **2))
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

        # prepare integration steps qOffset:
        qOffset = np.logspace(np.log10(q.min() / 10.),
                np.log10(xb / 2.), num = n)
        qOffset = np.concatenate(([0,], qOffset)) [np.newaxis, :]

        if xb == xt:
            y = 1. - (qOffset * 0.)
        else:
            y = 1. - (qOffset - xt) / (xb - xt)

        y = np.clip(y, 0., 1.)
        y[qOffset < xt] = 1.
        Area = (xt + 0.5 * (xb - xt))
        self._qOffset, self._weights = qOffset, (y / Area)
    
    def updateQUnit(self, newUnit):
        assert isinstance(newUnit, ScatteringVector)
        self.umbra.setUnit(newUnit)
        self.penumbra.setUnit(newUnit)

    def updateQLimits(self, qLimit):
        qLow, qHigh = qLimit
        self.umbra.setValueRange(qLimit)
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
        #print >>sys.__stderr__, "integrate", xb, xt

        # ensure things are what they are supposed to be
        assert (xt >= 0.)
        if xb < xt:
            xb = xt # should use square profile in this case.

        # prepare integration steps qOffset:
        qOffset = numpy.logspace(numpy.log10(q.min() / 10.),
                            numpy.log10(xb / 2.), num = n)
        qOffset = numpy.concatenate(([0,], qOffset)) [numpy.newaxis, :]

        if xb == xt: 
            y = 1. - (qOffset * 0.)
        else:
            y = 1. - (qOffset - xt) / (xb - xt)

        y = numpy.clip(y, 0., 1.)
        y[qOffset < xt] = 1.
        area = (xt + 0.5 * (xb - xt))
        self._qOffset, self._weights = qOffset, y / area

TrapezoidSmearing.factory()

class CallbackRegistry(object):
    _callbacks = None # registered callbacks on certain events

    @abstractproperty
    def callbackSlots(self):
        raise NotImplementedError

    def register(self, what, func):
        # check for the correct number of arguments of func as well?
        assert isCallable(func)
        self._assertPurpose(what)
        if self._callbacks is None: # lazy init
            self._callbacks = dict()
        if what not in self._callbacks:
            self._callbacks[what] = []
        self._callbacks[what].append(func)

    def callback(self, what, *args, **kwargs):
        self._assertPurpose(what)
        if self._callbacks is None:
            return
        funcLst = self._callbacks.get(what, [])
        for func in funcLst:
            if not isCallable(func):
                continue
            func(*args, **kwargs)

    def _assertPurpose(self, what):
        assert what in self.callbackSlots, (
            "'{}' not in predefined callback slots '{}'"
            .format(what, self.callbackSlots))

class DataConfig(AlgorithmBase, CallbackRegistry):
    parameters = (
        Parameter("xLow", 0., unit = NoUnit(),
            displayName = "lower {x} cut-off",
            valueRange = (0., numpy.inf), decimals = 1),
        Parameter("xHigh", numpy.inf, unit = NoUnit(),
            displayName = "upper {x} cut-off",
            valueRange = (0., numpy.inf), decimals = 1),
    )

    @property
    def callbackSlots(self):
        return set(("xlimits",))

    def __init__(self):
        super(DataConfig, self).__init__()
        self.xLow.setOnValueUpdate(self.updateXLimits)
        self.xHigh.setOnValueUpdate(self.updateXLimits)

    @mixedmethod
    def updateXLimits(self):
        if not self.xLow() <= self.xHigh():
            temp = self.xLow()
            self.xLow.setValue(self.xHigh())
            self.xHigh.setValue(temp)
        self.callback("xlimits", (self.xLow(), self.xHigh()))

class SASConfig(DataConfig):
    # TODO: implement callbacks for each value pair & unit
    # set units already for use in parameter definitions and UI
    _iUnit = NoUnit()
    _qUnit = NoUnit()
    _pUnit = NoUnit()
    _smearing = None
    shortName = "SAS data configuration"

    @property
    def callbackSlots(self):
        return super(SASConfig, self).callbackSlots | set(("qunit", "punit", "iunit"))

    # not used yet, need a signal to get the available q-range of data from filelist to datawidget
    def setQRange(self, qMin, qMax):
        """Sets available range of loaded data."""
        self.qLow.setValueRange((qMin, qMax))
        self.qHigh.setValueRange((qMin, qMax))
        if self.smearing is None:
            return
        self.smearing.updateQLimits((self.qLow(), self.qHigh()))

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
        if self.smearing is None:
            return
        self.smearing.integrate(q)
        qOffset, weights = self.smearing.prepared
        #print >>sys.__stderr__, "prepareSmearing"
        #print >>sys.__stderr__, unicode(self)
        # calculate the intensities at sqrt(q**2 + qOffset **2)
        return numpy.sqrt(numpy.add.outer(q**2, qOffset[0,:]**2))

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
        self.register("qunit", self.xLow.setUnit)
        self.register("qunit", self.xHigh.setUnit)
        if self.smearing is not None:
            self.register("qunit", self.smearing.updateQUnit)
            self.register("xlimits", self.smearing.updateQLimits)
        self.iUnit = ScatteringIntensity(u"(m sr)⁻¹")
        self.qUnit = ScatteringVector(u"nm⁻¹")
        self.pUnit = Angle(u"°")

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
