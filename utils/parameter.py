# -*- coding: utf-8 -*-
# utils/parameter.py

from cutesnake.utils import mixedmethod
from cutesnake.utils.tests import testfor, isInteger
from cutesnake.algorithm import (ParameterBase, ParameterFloat,
                                 ParameterNumerical, ParameterBoolean,
                                 ParameterLog, ParameterString)
from cutesnake.algorithm import Parameter as ParameterFactory

# put this sketch here, for the moment can be placed in a separate file later
class Histogram(object):
    """Stores histogram related settings of a parameter.
    The results too, eventually(?)
    """
    _binCount  = None
    _scaleX    = None
    _weighting = None

    @staticmethod
    def availableScaling(index = None):
        avail = ('lin', 'log')
        try:
            return avail[index]
        except:
            return avail

    @staticmethod
    def availableWeighting(index = None):
        avail = ('vol', 'num')
        try:
            return avail[index]
        except:
            return avail

    @property
    def binCount(self):
        return self._binCount

    @binCount.setter
    def binCount(self, binCount):
        testfor(isInteger(binCount) and binCount > 0,
                ValueError,
                "Histogram bin count has to be an integer larger 0!")
        self._binCount = max(0, int(binCount))

    @property
    def scaleX(self):
        return self._scaleX

    @scaleX.setter
    def scaleX(self, kind):
        self._scaleX = str(kind)
        if self._scaleX not in self.availableScaling():
            self._scaleX = self.availableScaling(0)

    @property
    def weighting(self):
        return self._weighting

    @weighting.setter
    def weighting(self, kind):
        self._weighting = str(kind)
        if self._weighting not in self.availableWeighting():
            self._weighting = self.availableWeighting(0)

    def __init__(self, binCount = 10):
        """Creates an histogram with default bin count, will be updated later."""
        self.binCount = binCount # bin count is mandatory
        self.scaleX = None # sets it to the first available option by default
        self.weighting = None # sets it to the first available option by def.

    # TODO: Function toString() or toJson() or similar which makes it
    # serializable and also a classmethod which constructs it from serial data

class FitParameterBase(ParameterBase):
    ParameterBase.setAttributes(locals(), isActive = False)



class FitParameterString(FitParameterBase, ParameterString):
    pass

class FitParameterBoolean(FitParameterBase, ParameterBoolean):
    pass

class FitParameterNumerical(FitParameterBase, ParameterNumerical):
    pass

class FitParameterFloat(FitParameterBase, ParameterFloat):
    pass

class FitParameterLog(FitParameterBase, ParameterLog):
    pass

def Parameter(*args, **kwargs):
    return ParameterFactory(*args, paramTypes = (
                FitParameterBoolean,
                FitParameterFloat,
                FitParameterNumerical,
                FitParameterString,
                FitParameterBase
            ), **kwargs)

# vim: set ts=4 sts=4 sw=4 tw=0:
