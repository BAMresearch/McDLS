# -*- coding: utf-8 -*-
# models/scatteringmodel.py

from abc import ABCMeta, abstractmethod
from cutesnake.utils import isList
from cutesnake.algorithm import AlgorithmBase
from utils.propertynames import PropertyNames

class ScatteringModel(AlgorithmBase, PropertyNames):
    __metaclass__ = ABCMeta
    compensationExponent = 0.5 # default

    def updateParamBounds(self, bounds):
        if not isList(bounds):
            bounds = [bounds,]
        if not isinstance(bounds, list):
            bounds = list(bounds)
        return bounds

    # it doesn't belong to the model?
    # should be instrumentation geometry ...
    def smear(self, arg):
        return arg

    @abstractmethod
    def vol(self, paramValues, compensationExponent = None):
        """Calculates the volume of this model, taking compensationExponent
        into account from input or preset parameters."""
        if len(self) == 0 and paramValues is None:
            return True
        return paramValues.shape[1] == len(self)

    @abstractmethod
    def ff(self, dataset, paramValues):
        """Calculates the Rayleigh function of this model."""
        if len(self) == 0 and paramValues is None:
            return True
        return paramValues.shape[1] == len(self)

    def generateParameters(self, count = 1):
        """Generates a set of parameters for this model using the predefined
        Parameter.generator. Allows for different random number distributions.
        """
        lst = numpy.zeros((count, len(self)))
        for idx, param in self:
            # generate numbers in different range for each parameter
            lst[:, idx] = param.generate(count = count)
        # output count-by-nParameters array
        return lst

# vim: set ts=4 sts=4 sw=4 tw=0:
