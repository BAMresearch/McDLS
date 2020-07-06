# -*- coding: utf-8 -*-
# bases/model/dlsmodeldata.py

from numpy import sqrt
from . import ModelData

class DLSModelData(ModelData):

    @property
    def chisqrInt(self):
        """Normalize and square the cumulated model intensities before passing
        them to the chi-square test. This is g1(tau)², the first-order
        correlation function squared."""
        return (self.cumInt / sum(self.wset))**2

    def volumeFraction(self, scaling):
        """Using the square root of the scaling factor to determine the volume
        fraction because the model intensities are squared after cumulation and
        normalization during post-processing."""
        return super(DLSModelData, self).volumeFraction(sqrt(scaling))

class DLSModelDataPlainVol(DLSModelData):
    def volumeFraction(self, scaling):
        return (super(DLSModelDataPlainVol, self)
                        .volumeFraction(scaling) * self.vset)

# vim: set ts=4 sts=4 sw=4 tw=0:
