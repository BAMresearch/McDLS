# -*- coding: utf-8 -*-
# dataobj/sasdata.py

"""
Represents data associated with a measurement by small angle scattering (SAS).
Some examples and tests.

>>> import numpy
>>> testdata = numpy.random.rand(4,4)
>>> testtitle = "some title"
>>> from sasdata import SASData

Testing
>>> first = SASData(testtitle, testdata)
>>> first.title == testtitle
True
>>> numpy.all(first.rawArray == testdata)
True
"""

import logging
import numpy as np # For arrays

from ..utils import classproperty
from ..utils.units import Length, ScatteringVector, ScatteringIntensity, Angle
from . import DataObj, SASConfig, DataVector
from ..bases.model import SASModel

class SASData(DataObj):
    """Represents one set of data from a unique source (a file, for example).
    """
    _e = None # internal DataVector
    _sizeEst = None
    _shannonChannelEst = None
    _rUnit = None # defines units for r used in sizeest

    # define DataObj interface

    @classproperty
    @classmethod
    def sourceName(cls):
        """The type of data source for UI label text."""
        return "Small Angle Scattering"

    @property
    def qLimsString(self):
        """Properly formatted q-limits for UI label text."""
        return self.x0.limsString

    @property
    def q(self):
        """Q-Vector at which the intensities are measured.
        Provided for convenience use within models."""
        return self.x0.binnedData # reverts to sanitized if not binned

    @property
    def pLimsString(self):
        """Properly formatted q-limits for UI label text."""
        if self.x1 is None:
            return ""
        return self.x1.limsString

    @property
    def p(self):
        """Q-Vector at which the intensities are measured.
        Provided for convenience use within models."""
        if self.x1 is None:
            return np.array(())
        return self.x1.sanitized

    # general information on this data set

    @property
    def count(self):
        return len(self.x0.binnedData) # used to be sanitized

    # general info texts for the UI

    @property
    def dataContent(self):
        """shows the content of the loaded data: Q, I, IErr, etc"""
        content = []
        if self.x0 is not None:
            content.append(self.x0.name)
        if self.f is not None:
            content.append(self.f.name)
        # if self.fu is not None: # cannot be none
        #     content.append(u"σI")
        if self.is2d:
            content.append('Psi')
        return ", ".join(content)

    @classproperty
    @classmethod
    def displayDataDescr(cls):
        return ("Filename ", "Data points ", "Data content ", 
                "Q limits ", "Est. sphere size ", "Recommended number of bins ")

    @classproperty
    @classmethod
    def displayData(cls):
        return ("title", "count", "dataContent", 
                "qLimsString", "sphericalSizeEstText", "shannonChannelEstText")

    def sphericalSizeEst(self):
        return self._sizeEst

    @property
    def rUnit(self):
        return self._rUnit

    @property
    def sphericalSizeEstText(self):
        return u"{0:.3g} ≤ R ({rUnitName}) ≤ {1:.3g} ".format(
                *self.rUnit.toDisplay(self.sphericalSizeEst()),
                rUnitName = self.rUnit.displayMagnitudeName)

    def shannonChannelEst(self):
        if self._shannonChannelEst is None:
            return None # otherwise int() below fails with TypeError
        return int(self._shannonChannelEst)

    def _prepareShannonChannelEst(self, *args):
        self._shannonChannelEst = self.x0.limit[1] / self.x0.limit[0]

    @property
    def shannonChannelEstText(self):
        est = self.shannonChannelEst()
        if est is None:
            return None
        return u"≤ {0:2d} bins ".format(est)

    def __init__(self, **kwargs):
        super(SASData, self).__init__(**kwargs)
        # self._h5LocAdd = "sasdata01" # overwriting DataObj default; DOES NOTHING
        
        # process rawArray for new DataVector instances:
        rawArray = kwargs.pop('rawArray', None)
        if rawArray is None:
            logging.error('SASData must be called with a rawArray provided')

        self.x0 = DataVector(u'q', rawArray[:, 0],
                             unit = ScatteringVector(u"nm⁻¹"))
        self.f  = DataVector(u'I', rawArray[:, 1], rawU = rawArray[:, 2],
                             unit = ScatteringIntensity(u"(m sr)⁻¹"))
        # sanitized uncertainty, we should use self._e.copy
        logging.info("Init SASData: " + self.qLimsString)
        if (rawArray.shape[1] > 3 and rawArray[:, 3].min()
                                   != rawArray[:, 3].max()):
            # psi column is present and contains some data
            self.x1 = DataVector(u'ψ', rawArray[:, 3], unit = Angle(u"°"))
            logging.info(self.pLimsString)

        #set unit definitions for display and internal units
        self._rUnit = Length(u"nm")
        # init config as early as possible to get properties ready which
        # depend on it (qlow/qhigh?)
        # (should be moved to DataObj but the DataVectors have to be set earlier)
        self.initConfig()

    def updateConfig(self):
        super(SASData, self).updateConfig()
        self.config.register("x0limits", self._prepareShannonChannelEst)
        # prepare
        self.locs = self.config.prepareSmearing(self.x0.binnedData)
        # suggested upgrade for 2d smearing:
        # self.locs = self.config.prepareSmearing(
        #                           self.x0.siData, self.x1.siData)

    @property
    def configType(self):
        return SASConfig

    @property
    def modelType(self):
        return SASModel

    def _propagateMask(self, *args):
        super(SASData, self)._propagateMask(*args)
        if self.x0.limit[0] == 0.:
            return
        self._sizeEst = np.pi / np.array([self.x0.limit[1],
                                          abs(self.x0.limit[0])])

if __name__ == "__main__":
    import doctest
    doctest.testmod()

# vim: set ts=4 sts=4 sw=4 tw=0:
