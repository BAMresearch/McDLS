# -*- coding: utf-8 -*-
# bases/model/dlsmodel.py

from abc import ABCMeta
from bases.algorithm import Parameter
from utils.units import NoUnit
from bases.model import ScatteringModel, DLSModelData, DLSModelDataPlainVol

class DLSModel(ScatteringModel):
    __metaclass__ = ABCMeta
    _scatteringVector = None
    parameters = (
        Parameter("ampSquared", True, unit = NoUnit(),
            displayName = "square the amplitude?",
            description = "ON = (volume * formfactor)², "
                          "OFF = volume * formfactor"),
    )

    def modelDataType(self):
        if self.ampSquared():
            return DLSModelData
        else:
            return DLSModelDataPlainVol

    @property
    def scatteringVector(self):
        return self._scatteringVector

    def calcIntensity(self, data, **kwargs):
        self._scatteringVector = data.scatteringVector
        v = self._volume()
        w = self._weight()
        s = self.surface()
        if self.ampSquared():
            w *= w # square the weight, i.e. amplitude
        # calculate their form factors
        ff = self._formfactor(data)
        # a set of intensities
        it = ff * w
        return it, v, w, s

# vim: set ts=4 sts=4 sw=4 tw=0:
