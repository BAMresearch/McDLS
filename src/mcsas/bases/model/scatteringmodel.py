# -*- coding: utf-8 -*-
# bases/model/scatteringmodel.py

import os.path
import logging
from abc import ABCMeta, abstractmethod
from future.utils import with_metaclass
from numpy import arange, zeros, argmax, hstack

from ...utils import isList, mixedmethod, testfor, classname
from ..algorithm import AlgorithmBase
from ...utils.parameter import isActiveFitParam

class ScatteringModel(with_metaclass(ABCMeta, AlgorithmBase)):
    @abstractmethod
    def volume(self):
        """Calculates the volume of this model, taking compensationExponent
        into account from input or preset parameters.
        Reimplement this for new models."""
        raise NotImplementedError

    def absVolume(self):
        """Forwarding to usual volume() by default.
        Can be overridden to include SLD."""
        return self.volume()

    def _volume(self, compensationExponent = None):
        """Wrapper around the user-defined function."""
        self.compensationExponent = compensationExponent
        # calling user provided custom model
        return self.absVolume()

    @abstractmethod
    def weight(self):
        """A weighting function for the form factor.
        With SAXS, it is usually the volume squared."""
        raise NotImplementedError

    def _weight(self, compensationExponent = None):
        """Wrapper around the user-defined function."""
        self.compensationExponent = compensationExponent
        return self.weight()

    @abstractmethod
    def formfactor(self, dataset):
        """Calculates the Rayleigh function of this model.
        Reimplement this for new models."""
        raise NotImplementedError

    def _formfactor(self, dataset):
        """Wrapper around the user-defined function."""
        # calling user provided custom model
        i = self.formfactor(dataset)
        return i

    def surface(self):
        """Returns the surface area of a single scatterer. Used for the
        surface weighted distribution histogram. Returns 0 by default.
        Reimplement this for a model."""
        return 0

    def hdfWrite(self, hdf):
        hdf.writeAttributes(modelName = classname(self))
        for p in self.params():
            logging.debug("Writing model parameter: {} to HDF5".format(p.name()))
            hdf.writeMember(self, p.name())

    @abstractmethod
    def calcIntensity(self, data, compensationExponent = None):
        """Calculates the model intensity which is later compared to the data.
        Returns a tuple containing an array of the calculated intensities for
        the grid provided with the data and the volume of a single particle
        based on the model parameters.
        Has to be implemented in derived classes specific to a certain type of
        measurement.
        """
        raise NotImplementedError

    def calc(self, data, pset, compensationExponent = None):
        """Calculates the total intensity and scatterer volume contributions
        using the current model.
        *pset* number columns equals the number of active parameters.
        Returns a ModelData object for a certain type of measurement.
        """
        # remember parameter values
        params = self.activeParams()
        oldValues = [p() for p in params] # this sucks. But we dont want to lose the user provided value
        cumInt = zeros(data.f.binnedData.shape) # cumulated intensities
        vset = zeros(pset.shape[0])
        wset = zeros(pset.shape[0])
        sset = zeros(pset.shape[0])
        # call the model for each parameter set explicitly
        # otherwise the model gets complex for multiple params incl. fitting
        for i in arange(pset.shape[0]): # for each contribution
            for p, v in zip(params, pset[i]): # for each fit param within
                p.setValue(v)
            # result squared or not is model type dependent
            it, vset[i], wset[i], sset[i] = self.calcIntensity(data,
                          compensationExponent = compensationExponent)
            # a set of intensities
            cumInt += it
        # restore previous parameter values
        for p, v in zip(params, oldValues):
            p.setValue(v)
        return self.getModelData(cumInt, vset, wset, sset)

    def getModelData(self, cumInt, vset, wset, sset):
        return self.modelDataType()(cumInt.flatten(), vset, wset, sset,
                                    self.activeParamCount())

    @abstractmethod
    def modelDataType(self):
        """Returns the appropriate ModelData class for this type of model.
        """
        raise NotImplementedError

    def generateParameters(self, count = 1):
        """Generates a set of parameters for this model using the predefined
        Parameter.generator. Allows for different random number distributions.
        """
        lst = zeros((count, self.activeParamCount()))
        for idx, param in enumerate(self.activeParams()):
            # generate numbers in different range for each active parameter
            if isActiveFitParam(param):
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
        aPars = [par for par in setforcls.params() if isActiveFitParam(par)]
        return tuple(aPars)

    @mixedmethod
    def activeParamCount(setforcls):
        return len(setforcls.activeParams())

    @mixedmethod
    def activeParamNames(setforcls):
        namelist = list()
        for param in setforcls.activeParams():
            namelist.append(param.displayName())
        return namelist

    # helpers for model testing below

    @mixedmethod
    def update(selforcls, paramDict):
        """Update parameter values based on provided dict with parameter
        names as keys."""
        selforcls.fixTestParams(paramDict)
        for key, value in paramDict.items():
            try:
                p = getattr(selforcls, key)
            except: pass
            else:
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
        # exclude the base name at front
        pnames = tuple(enumerate(pnames[1:]))
        # store all values by index for reference in fixTestParams()
        result = dict(pnames)
        # map values to parameter names beginning at front
        for i, valueStr in pnames:
            try:
                p = cls.param(i) # raises IndexError eventually
                result[p.name()] = p.dtype(valueStr)
            except IndexError: continue
        return result

    @classmethod
    def test(cls, filename):
        """Regression test of a scattering model. File names are expected
        to contain the parameter values which produce the provided intensity.
        Otherwise implement fixTestParams() for the particular model.

        - *filename*: Name of the file in cls.testDataDir to test against.
        - *cls.testRelErr*: Acceptable mean of relative error against reference
                            intensity. Default: 1e-5
        - *cls.testVolExp*: Volume compensation exponent, sets the amount of
                            volume contribution the intensity is scaled by.
        - *cls.testDataDir*: Directory of test data relative to program dir.
                             Default: "testdata"
        """
        return # this does not work anymore
        relerr = getattr(cls, "testRelErr", 1e-5)
        datadir = getattr(cls, "testDataDir", "testdata")
        volumeExponent = getattr(cls, "testVolExp", 1.0)
        filename = os.path.abspath(os.path.join(datadir, filename))
        #dataset = SASData.load(filename)
        if dataset is None:
            return
        model = cls()
        testParams = model.getParametersFromFilename(dataset.filename)
        model.update(testParams)
        # intensity how it's calculated in SASfit
        intensity = (model.vol(None,
                               compensationExponent = volumeExponent)
                     * model._formfactor(dataset, None))**2.
        # computing the relative error to reference data
        delta = abs((dataset.f.binnedData - intensity) / dataset.f.binnedData)
        dmax = argmax(delta)
        testfor(delta.mean() < relerr, AssertionError,
                "Could not verify {model} intensity against\n'{fn}',"
                "\nmean(delta): {mean} >= relErr: {relerr}"
                "\nQ, Int_ref, Int_calc, Delta around dmax({dmax}):"
                "\n{data}"
                .format(model = cls.name(), fn = filename,
                        mean = delta.mean(), relerr = relerr,
                        dmax = dmax, data = hstack((
                            dataset.x0.binnedData.reshape(-1, 1),
                            dataset.f.binnedData.reshape(-1, 1),
                            intensity.reshape(-1, 1),
                            delta.reshape(-1, 1)))[max(0, dmax-4):dmax+5]
                        )
                )

# vim: set ts=4 sts=4 sw=4 tw=0:
