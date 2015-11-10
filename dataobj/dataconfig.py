# -*- coding: utf-8 -*-
# dataobj/dataconfig.py

from __future__ import absolute_import # PEP328

from abc import ABCMeta, abstractmethod
import numpy
from bases.algorithm import AlgorithmBase
from utils.parameter import Parameter
from utils.mixedmethod import mixedmethod
from utils.units import ScatteringIntensity, ScatteringVector, Angle, NoUnit
from utils import clip

class SmearingConfig(AlgorithmBase):
    """Abstract base class, can't be instantiated."""
    __metaclass__ = ABCMeta
    parameters = (
        # not sure if this is the right place: is the nsteps parameter useful
        # for all possible smearing settings?
        Parameter("nSteps", 25, unit = NoUnit(),
            displayName = "number of smearing points around each q",
            valueRange = (0, 1000)),
    )

    @abstractmethod
    def updateQUnit(self, newUnit):
        pass

    @abstractmethod
    def updateQLimits(self, qLow, qHigh):
        pass

import sys

class TrapezoidSmearing(SmearingConfig):
    parameters = (
        Parameter("umbra", 0., unit = NoUnit(), # unit set outside
            displayName = "top width of the trapezoidal beam length profile",
            valueRange = (0., numpy.inf), decimals = 1),
        Parameter("penumbra", 0., unit = NoUnit(), # unit set outside
            displayName = "bottom width of the trapezoidal beam length profile",
            valueRange = (0., numpy.inf), decimals = 1),
    )

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
        print >>sys.__stderr__, "onUmbraUpdate"
        # value in Penumbra must not be smaller than Umbra
        self.penumbra.setValueRange((self.umbra(), self.penumbra.max()))

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

    def __init__(self):
        super(SASConfig, self).__init__()
        self.smearing = TrapezoidSmearing()
        self.iUnit = ScatteringIntensity(u"(m sr)⁻¹")
        self.qUnit = ScatteringVector(u"nm⁻¹")
        self.pUnit = Angle(u"°")
        self.qLow.setOnValueUpdate(self.updateConstraints)
        self.qHigh.setOnValueUpdate(self.updateConstraints)

    def __str__(self):
        return "\n".join((
            super(SASConfig, self).__str__(),
            unicode(self.smearing)
        ))

SASConfig.factory() # check class variables

# vim: set ts=4 sts=4 sw=4 tw=0:
