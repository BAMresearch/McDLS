# -*- coding: utf-8 -*-
# dataobj/dataconfig.py

from __future__ import absolute_import # PEP328

from abc import ABCMeta, abstractproperty
import numpy
from bases.algorithm import AlgorithmBase
from utils.parameter import Parameter
from utils.mixedmethod import mixedmethod
from utils.units import NoUnit
from utils import isCallable

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
        if func not in self._callbacks[what]:
            self._callbacks[what].append(func)

    def callback(self, what, *args, **kwargs):
        self._assertPurpose(what)
        if self._callbacks is None:
            return
        funcLst = []
        for func in self._callbacks.get(what, []):
            if not isCallable(func):
                continue
            func(*args, **kwargs)
            funcLst.append(func)
        # update the callback list, invalid functions removed
        self._callbacks[what] = funcLst

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

# vim: set ts=4 sts=4 sw=4 tw=0:
