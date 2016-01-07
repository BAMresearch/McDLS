# -*- coding: utf-8 -*-
# dataobj/dataconfig.py

from __future__ import absolute_import # PEP328

from abc import ABCMeta, abstractproperty
from types import MethodType
import numpy
from bases.algorithm import AlgorithmBase
from utils.parameter import Parameter
from utils.mixedmethod import mixedmethod
from utils.units import NoUnit
from utils import isCallable

def funcNotInFuncList(f, flst):
    """Custom predicate for comparing bounded methods:
    Duplicate only if instance ID and method name match.
    """
    if not isinstance(f, MethodType):
        return (f not in flst)

    idExists = (id(f.im_self) in [id(of.im_self)
                    for of in flst if isinstance(of, MethodType)])
    if not idExists:
        return True
    nameExists = (f.im_func.func_name in [of.im_func.func_name
                    for of in flst if isinstance(of, MethodType)])
    if not nameExists:
        return True
    return False # both exist

class CallbackRegistry(object):
    _callbacks = None # registered callbacks on certain events

    @abstractproperty
    def callbackSlots(self):
        raise NotImplementedError

    def register(self, what, *func):
        # check for the correct number of arguments of func as well?
        assert all((isCallable(f) for f in func))
        self._assertPurpose(what)
        if self._callbacks is None: # lazy init
            self._callbacks = dict()
        if what not in self._callbacks:
            self._callbacks[what] = []
        for f in func:
            if funcNotInFuncList(f, self._callbacks[what]):
                self._callbacks[what].append(f)

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

    def __getstate__(self):
        state = self.__dict__.copy()
        state["_callbacks"] = None
        return state

class DataConfig(AlgorithmBase, CallbackRegistry):
    _is2d = False
    parameters = (
        Parameter("x0Low", 0., unit = NoUnit(),
            displayName = "lower {x0} cut-off",
            valueRange = (0., numpy.inf), decimals = 1),
        Parameter("x0LowClip", 0, unit = NoUnit(),
            displayName = "ignore leading {x0} points",
            valueRange = (0, 2**30)),
        Parameter("x0High", numpy.inf, unit = NoUnit(),
            displayName = "upper {x0} cut-off",
            valueRange = (0., numpy.inf), decimals = 1),
        Parameter("x1Low", 0., unit = NoUnit(),
            displayName = "lower {x1} cut-off",
            valueRange = (0., numpy.inf), decimals = 1),
        Parameter("x1High", numpy.inf, unit = NoUnit(),
            displayName = "upper {x1} cut-off",
            valueRange = (0., numpy.inf), decimals = 1),
        Parameter("fMaskZero", False, unit = NoUnit(),
            displayName = "Mask {f} values of 0", description =
            "Renders intensity values that are zero invalid for fitting"),
        Parameter("fMaskNeg", False, unit = NoUnit(),
            displayName = "Mask negative {f} values", description =
            "Renders negative intensity values invalid for fitting"),
    )

    @property
    def showParams(self):
        lst = super(DataConfig, self).showParams
        if not self.is2d: # put psi settings right behind q settings
            lst.remove("x1Low")
            lst.remove("x1High")
        return lst

    @property
    def callbackSlots(self):
        return set(("x0limits", "x1limits", "x0Clipping", "fMasks"))

    @property
    def is2d(self):
        return self._is2d

    @is2d.setter
    def is2d(self, isit):
        self._is2d = isit

    def __init__(self):
        super(DataConfig, self).__init__()
        self.x0Low.setOnValueUpdate(self.updateX0Limits)
        self.x0LowClip.setOnValueUpdate(self.updateX0Clipping)
        self.x0High.setOnValueUpdate(self.updateX0Limits)
        self.x1Low.setOnValueUpdate(self.updateX1Limits)
        self.x1High.setOnValueUpdate(self.updateX1Limits)
        self.fMaskZero.setOnValueUpdate(self.updateFMasks)
        self.fMaskNeg.setOnValueUpdate(self.updateFMasks)

    def updateX0Limits(self):
        self._onLimitUpdate("x0limits", self.x0Low, self.x0High)

    def updateX1Limits(self):
        self._onLimitUpdate("x1limits", self.x1Low, self.x1High)

    def _onLimitUpdate(self, callbackName, pLow, pHigh):
        if not pLow() <= pHigh():
            temp = pLow()
            pLow.setValue(pHigh())
            pHigh.setValue(temp)
        self.callback(callbackName, (pLow(), pHigh()))

    def updateX0Clipping(self):
        self.callback("x0Clipping", self.x0LowClip())

    def updateFMasks(self):
        self.callback("fMasks", (self.fMaskZero(), self.fMaskNeg()))

    def setX0ValueRange(self, limits):
        """Sets available range of loaded data."""
        self.x0Low.setValueRange(limits)
        self.x0High.setValueRange(limits)

    def setX1ValueRange(self, limits):
        pass

    def __getstate__(self):
        view = self.__dict__.viewkeys()
        state = dict()
        for cls in type(self).mro():
            if issubclass(cls, DataConfig) or not hasattr(cls, "__getstate__"):
                continue
            parentState = cls.__getstate__(self)
            view = view & parentState.viewkeys()
            state.update([(key, parentState[key]) for key in view])
        return state

    def __setstate__(self, state):
        super(DataConfig, self).__setstate__(state)

DataConfig.factory()

# vim: set ts=4 sts=4 sw=4 tw=0:
