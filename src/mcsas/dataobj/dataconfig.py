# -*- coding: utf-8 -*-
# dataobj/dataconfig.py


from builtins import object
from abc import ABCMeta, abstractproperty
from types import MethodType
import numpy
from bases.algorithm import AlgorithmBase
from bases.algorithm import Parameter # not defined in utils.parameter
from utils.mixedmethod import mixedmethod
from utils.units import NoUnit, Fraction
from utils import isCallable, classname

def funcNotInFuncList(f, flst):
    """Custom predicate for comparing bounded methods:
    Duplicate only if instance ID and method name match.
    """
    if not isinstance(f, MethodType):
        return (f not in flst)

    idExists = (id(f.__self__) in [id(of.__self__)
                    for of in flst if isinstance(of, MethodType)])
    if not idExists:
        return True
    nameExists = (f.__func__.__name__ in [of.__func__.__name__
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
    _sampleName = None
    _x0seen, _x1seen = None, None # remembers data sets seen
    parameters = (
        Parameter("x0Low", 0., unit = NoUnit(),
            displayName = "lower {x0} cut-off",
            valueRange = (0., numpy.inf), decimals = 10),
        Parameter("x0High", numpy.inf, unit = NoUnit(),
            displayName = "upper {x0} cut-off",
            valueRange = (0., numpy.inf), decimals = 10),
        Parameter("x1Low", 0., unit = NoUnit(),
            displayName = "lower {x1} cut-off",
            valueRange = (0., numpy.inf), decimals = 10),
        Parameter("x1High", numpy.inf, unit = NoUnit(),
            displayName = "upper {x1} cut-off",
            valueRange = (0., numpy.inf), decimals = 10),
        Parameter("fMaskZero", False, unit = NoUnit(),
            displayName = "Mask {f} values of 0", description =
            "Renders intensity values that are zero invalid for fitting"),
        Parameter("fMaskNeg", False, unit = NoUnit(),
            displayName = "Mask negative {f} values", description =
            "Renders negative intensity values invalid for fitting"),
        Parameter("fuMin", Fraction(u"%").toSi(1.), unit = Fraction(u"%"),
            displayName = "minimum uncertainty estimate", description =
            "Ensures that the uncertainties of the data are greater than the "
            "given fraction of the measured signal.",
            valueRange = (0., 1.), decimals = 9),
        Parameter("fuConst", False, unit = NoUnit(),
            displayName = "set uncertainties constant", description =
            "Sets all uncertainties of the data to the same value, "
            "such as 0.1, resulting in equal weighting of all data points."
            "(mutually exclusive with relative uncertainty)"),
        Parameter("fuRel", False, unit = NoUnit(),
            displayName = u"use relative uncertainty [σ/{f}]", description =
            "Divides all uncertainties of the data by the corresponding "
            "data value to weigh it relatively."
            "(mutually exclusive with constant uncertainty)"),
        Parameter("nBin", 100, unit = NoUnit(),
            displayName = "target number of bins \n (0 = no re-binning)",
            description = "Sets the number of bins to rebin the data into. \n"
                          "May be smaller than target value.",
            valueRange = (0., 1000)),
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
        return set(("x0limits", "x1limits", "fMasks", "fuMin"))

    @property
    def is2d(self):
        return self._is2d

    @is2d.setter
    def is2d(self, isit):
        self._is2d = isit

    @property
    def sampleName(self):
        return self._sampleName

    @sampleName.setter
    def sampleName(self, newName):
        self._sampleName = newName

    def __init__(self):
        super(DataConfig, self).__init__()
        self.x0Low.setOnValueUpdate(self.updateX0Limits)
        self.x0High.setOnValueUpdate(self.updateX0Limits)
        self.x1Low.setOnValueUpdate(self.updateX1Limits)
        self.x1High.setOnValueUpdate(self.updateX1Limits)
        self.fMaskZero.setOnValueUpdate(self.updateFMasks)
        self.fMaskNeg.setOnValueUpdate(self.updateFMasks)
        self.fuMin.setOnValueUpdate(self.updateFuMin)
        self.fuConst.setOnValueUpdate(self.updateFuConst)
        self.fuRel.setOnValueUpdate(self.updateFuRel)

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

    def updateFMasks(self):
        self.callback("fMasks", (self.fMaskZero(), self.fMaskNeg()))

    def updateFuMin(self):
        self.callback("fuMin", self.fuMin())

    def updateFuConst(self):
        if self.fuConst():
            self.fuRel.setValue(False)
        # disable minimum uncertainty to avoid accidental overwriting
        self.fuMin.setValue(0.0)
        # hand over to fuMin
        self.updateFuMin()

    def updateFuRel(self):
        if self.fuRel():
            self.fuConst.setValue(False)
        # disable minimum uncertainty to avoid accidental overwriting
        self.fuMin.setValue(0.0)
        # hand over to fuMin
        self.updateFuMin()

    def updateX0Unit(self, newUnit):
        """Sets the unit of the x0 vector."""
        # No callback this time because it is updated top-down from the
        # DataVectors in DataObj. The unit is not expected the be updated
        # in the Parameter only and has to bubble upwards to the DataVector
        # (atm!).
        self.x0Low.setUnit(newUnit)
        self.x0High.setUnit(newUnit)

    def updateX1Unit(self, newUnit):
        self.x1Low.setUnit(newUnit)
        self.x1High.setUnit(newUnit)

    def onUpdatedX0(self, x0):
        """Sets available range of loaded data."""
        if self._x0seen is None:
            # on the first data, set the param limits to the exact value range
            limits = (x0.min(), x0.max())
            self._x0seen = id(x0) # just store something for now
        else: # there were other data sets already, the value range grows
            # alternatives: (1) shrinking means another range for broader
            # datasets can not be selected in the UI;
            limits = self.x0Low.valueRange()
            limits = min(x0.min(), limits[0]), max(x0.max(), limits[1])
        self.x0Low.setValueRange(limits)
        self.x0High.setValueRange(limits)

    def onUpdatedX1(self, x1):
        pass

    def hdfWrite(self, hdf):
        super(DataConfig, self).hdfWrite(hdf)
        hdf.writeMembers(self, 'sampleName', 'is2d')

# vim: set ts=4 sts=4 sw=4 tw=0:
