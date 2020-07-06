# -*- coding: utf-8 -*-
# bases/model/sasmodel.py

import logging
import numpy as np
from abc import ABCMeta
from future.utils import with_metaclass

from . import SASModelData, ScatteringModel

class SASModel(with_metaclass(ABCMeta, ScatteringModel)):
    canSmear = False # Indicates a model function which supports smearing...

    def modelDataType(self):
        return SASModelData

    def __init__(self):
        # just checking:
        super(SASModel, self).__init__()
        logging.debug("SASData init method called")

    def getQ(self, dataset):
        """ This is a function that returns Q. In case of smearing, dataset itself
        is a 2D matrix of Q-values. When smearing is not enabled, dataset.q contains
        a 1D vector of q.

        I do realize that this is not a good way of doing things. This should be
        replaced at a given point in time by a better solution within sasdata.
        """

        if isinstance(dataset, np.ndarray):
            q = dataset
        else:
            q = dataset.q
        return q

    def weight(self):
        r"""Calculates an intensity weighting used during fitting. It is based
        on the scatterers volume. It can be modified by a user-defined
        compensation exponent *c*. The default value is :math:`c={2 \over 3}`

        :math:`w(r) = v(r)^{2c}`
        """
        return self.volume()**(2 * self.compensationExponent)

    def calcIntensity(self, data, compensationExponent = None):
        r"""Returns the intensity *I*, the volume :math:`v_{abs}` and the
        intensity weights *w* for a single parameter contribution over all *q*:

        :math:`I(q,r) = F^2(q,r) \cdot w(r)`
        """
        v = self._volume(compensationExponent = compensationExponent)
        w = self._weight(compensationExponent = compensationExponent)
        s = self.surface()

        if ((data.config.smearing is not None) and
                self.canSmear and
                data.config.smearing.doSmear() and # serves same purpose as first
                data.config.smearing.inputValid()):
            # inputValid can be removed once more appropriate limits are set in GUI

            # TODO: fix after change from x0Fit to x0:
            locs = data.locs # [data.x0.validIndices] # apply xlimits
            # the ff functions might only accept one-dimensional q arrays
            # kansas = locs.shape
            # locs = locs.reshape((locs.size))
            ff = self._formfactor(locs) # .reshape(kansas)
            qOffset, weightFunc = data.config.smearing.prepared
#            import sys
#            print >>sys.__stderr__, "prepared"
#            print >>sys.__stderr__, unicode(data.config.smearing)
            it = 2 * np.trapz(ff**2 * w * weightFunc,
                    x = qOffset, axis = 1)
        else:
            # calculate their form factors
            ff = self._formfactor(data)
            # a set of intensities
            it = ff**2 * w
        return it, v, w, s

# vim: set ts=4 sts=4 sw=4 tw=0:
