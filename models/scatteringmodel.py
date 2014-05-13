# -*- coding: utf-8 -*-
# models/scatteringmodel.py

import os.path
import logging
import inspect
from abc import ABCMeta, abstractmethod
import numpy
from utils import isList, isNumber, mixedmethod, testfor
from cutesnake.algorithm import AlgorithmBase
from utils.propertynames import PropertyNames
from sasdata import SASData

class ScatteringModel(AlgorithmBase, PropertyNames):
    __metaclass__ = ABCMeta
    compensationExponent = 0.5 # default

    # it doesn't belong to the model?
    # should be instrumentation geometry ...
    def smear(self, arg):
        return arg

    @abstractmethod
    def volume(self, paramValues):
        """Calculates the volume of this model, taking compensationExponent
        into account from input or preset parameters.
        Reimplement this for new models."""
        raise NotImplemented

    def vol(self, paramValues, compensationExponent = None):
        if paramValues is None:
            assert self.paramCount() == 0
        # by definition, no vector input, would require complex model code
        # compare the flat length
        assert paramValues.size == self.paramCount()
        self.compensationExponent = compensationExponent
        if self.compensationExponent is None:
            # get the exponent from class level of the particular model
            self.compensationExponent = type(self).compensationExponent
        # calling user provided custom model
        v = self.volume(paramValues)
        # volume always returns a single value
        assert isNumber(v)
        return v

    @abstractmethod
    def formfactor(self, dataset, paramValues = None):
        """Calculates the Rayleigh function of this model.
        Reimplement this for new models."""
        raise NotImplemented

    def ff(self, dataset, paramValues = None):
        if paramValues is None:
            assert self.paramCount() == 0
        # by definition, no vector input, would require complex model code
        # compare the flat length
        assert paramValues.size == self.paramCount()
        # calling user provided custom model
        i = self.formfactor(dataset, paramValues)
        # there has to be one intensity value for each q-vector
        assert i.size == dataset.q.size
        return i

    def generateParameters(self, count = 1):
        """Generates a set of parameters for this model using the predefined
        Parameter.generator. Allows for different random number distributions.
        """
        lst = numpy.zeros((count, self.paramCount()))
        for idx, param in enumerate(self.params()):
            # generate numbers in different range for each parameter
            #only for active parameters, otherwise it may try to generate
            #random values for a boolean-type parameter.
            if param.isActive():
                lst[:, idx] = param.generate(count = count)
        # output count-by-nParameters array
        return lst

    def updateParamBounds(self, bounds):
        if not isList(bounds):
            bounds = [bounds,]
        if not isinstance(bounds, list):
            bounds = list(bounds)
        return bounds

    @mixedmethod
    def activeParams(setforcls):
        """returns all "active" parameters of this algorithm"""
        aPars = [par for par in setforcls.params() if par.isActive()]
        return tuple(aPars)

    @mixedmethod
    def activeParamCount(setforcls):
        return len(setforcls.activeParams())

    @mixedmethod
    def update(selforcls, **kwargs):
        """Update parameter values based on provided dict with parameter
        names as keys."""
        selforcls.fixTestParams(kwargs)
        for key, value in kwargs.iteritems():
            p = getattr(selforcls, key, None)
            if p is None:
                continue
            p.setValue(p.dtype(value))

    @mixedmethod
    def fixTestParams(selforcls, params):
        """Eventually convert test parameters to the current model.
        Override this method in your model if parameters differ in order or
        meaning from those provided with the test data.
        *params*: A dict with mapping of parameter name to value.
        """
        return params

    @classmethod
    def getParametersFromFilename(cls, filename):
        """Derives model parameters for testing from reference data file."""
        pnames = os.path.splitext(os.path.basename(filename))[0]
        pnames = pnames.split("-")
        errorMsg = ("Could not infer {model} parameters from '{fn}'!"
                    .format(model = cls.name(), fn = filename))
        testfor(len(pnames) > cls.paramCount(), NameError, errorMsg)
        result = dict()
        try:
            for i, p in enumerate(cls.params()):
                value = pnames[i - cls.paramCount()]
                result[p.name()] = p.dtype(value)
        except:
            logging.error(errorMsg)
            raise
        return result

    @classmethod
    def test(cls, filenames = None, datadir = None, tol = 1e-5):
        """Regression test of a scattering model. File names are expected
        to contain the parameter values which produce the provided intensity.
        Otherwise implement fixTestParams() for the particular model.
        """
        testfor(filenames is not None and len(filenames), AssertionError,
                "No filenames for testing provided!")
        if datadir is None or not os.path.isdir(datadir):
            datadir = "../brianpauw"
        models = getModels()
        testfor(len(models), AssertionError,
                "No models found for testing!")
        for fn in filenames:
            fn = os.path.abspath(os.path.join(datadir, fn))
            dataset = SASData.load(fn)
            model = models[0]()
            testParams = model.getParametersFromFilename(dataset.filename)
            model.update(**testParams)
            intensity = model.formfactor(dataset, None)**2.
            # computing the relative error to reference data
            delta = abs((dataset.i - intensity) / dataset.i)
            testfor(delta.mean() < tol, AssertionError,
                    "Could not verify {model} intensity against\n'{fn}',\n"
                    "mean: {mean} >= tol: {tol}\n{delta}"
                    .format(model = cls.name(), fn = fn,
                            mean = delta.mean(), tol = tol, delta = delta))

def getModels():
    """Returns all subclasses of ScatteringModel in the current model script.
    Make sure classes are imported from top level modules in application path!
    """
    _dict = inspect.currentframe().f_back.f_back.f_globals
    return [obj for key, obj in _dict.items()
                if (isinstance(obj, type) and
                    issubclass(obj, ScatteringModel) and
                    obj is not ScatteringModel)]

# vim: set ts=4 sts=4 sw=4 tw=0:
