# -*- coding: utf-8 -*-
# mcsas/sasdata.py

"""
Represents data associated with a measurement by small angle scattering (SAS).
Some examples and tests.

>>> import numpy
>>> testdata = numpy.random.rand(4,4)
>>> testtitle = "some title"
>>> from sasdata import SASData

Testing copy()
>>> first = SASData(testtitle, testdata)
>>> first.title == testtitle
True
>>> numpy.all(first.origin == testdata)
True
>>> second = first.copy()
>>> second.title == testtitle
True
>>> numpy.all(second.origin == testdata)
True
>>> first == second
True
"""

from __future__ import absolute_import # PEP328
import os # Miscellaneous operating system interfaces
import logging
import numpy as np # For arrays
from bases.datafile import PDHFile, AsciiFile
from bases.dataset import DataSet, DisplayMixin
from utils import isList, classproperty
from gui.utils import processEventLoop
from utils.units import Length, ScatteringVector, ScatteringIntensity, Angle

class SASData(DataSet, DisplayMixin):
    """Represents one set of data from a unique source (a file, for example).
    """
    _sizeEst = None
    _eMin = 0.01 # minimum possible error (1%)
    _uncertainty = None
    _filename = None
    _qUnit = None # will be instance of SASUnit, defines units
    _pUnit = None # will be instance of SASUnit, defines units
    _iUnit = None # will be instance of SASUnit, defines units
    _rUnit = None # defines units for r used in sizeest
    _qClipRange = [-np.inf, np.inf] # Q clip range for this dataset
    _pClipRange = [-np.inf, np.inf] # psi clip range for this dataset
    _maskZeroInt = False # mask zero intensity values (updated by mcsasparam)
    _maskNegativeInt = False # mask negative intensity values (updated by mcsasparam)

    @classproperty
    @classmethod
    def displayDataDescr(cls):
        return ("filename", "data points", "data content", 
                "Q limits", "est. sphere size")

    @classproperty
    @classmethod
    def displayData(cls):
        return ("title", "count", "dataContent", 
                "qLimsString", "sphericalSizeEstText")

    @property
    def maskZeroInt(self):
        return bool(self._maskZeroInt)

    @maskZeroInt.setter
    def maskZeroInt(self, value):
        self._maskZeroInt = bool(value)

    @property
    def maskNegativeInt(self):
        return bool(self._maskNegativeInt)

    @maskNegativeInt.setter
    def maskNegativeInt(self, value):
        self._maskNegativeInt = bool(value)

    @property
    def count(self):
        return len(self.q)

    @property
    def sphericalSizeEstText(self):
        return u"{0:.3g} ≤ R ({rUnitName}) ≤ {1:.3g}".format(
                *self.rUnit.toDisplay(self.sphericalSizeEst()),
                rUnitName = self.rUnit.displayMagnitudeName)

    @property
    def eMin(self):
        return self._eMin

    @eMin.setter
    def eMin(self, value):
        value = float(value) 
        assert((value > 0.) and (value < 1.))
        self._eMin = value

    @property
    def qClipRange(self):
        return self._qClipRange

    @property
    def qMin(self):
        #returns minimum q from data or qClipRange, whichever is larger
        return np.maximum(self.qClipRange[0], self.q.min())

    @qMin.setter
    def qMin(self, newParam):
        #value in cliprange will not exceed available q.
        self.qClipRange[0] = np.maximum(newParam, self.q.min())

    @property
    def qMax(self):
        return np.minimum(self.qClipRange[1], self.q.max())

    @qMax.setter
    def qMax(self, newParam):
        #value in cliprange will not exceed available q.
        self.qClipRange[1] = np.minimum(newParam, self.q.max())

    @property
    def qLimsString(self):
        return u"{0:.3g} ≤ Q ({qMagnitudeName}) ≤ {1:.3g}".format(
                self.qUnit.toDisplay(self.qMin),
                self.qUnit.toDisplay(self.qMax),
                qMagnitudeName = self.qUnit.displayMagnitudeName)

    @property
    def pClipRange(self):
        return self._pClipRange

    @property
    def pMin(self):
        #returns minimum psi from data or psiClipRange, whichever is larger
        return np.maximum(self.pClipRange[0], self.p.min())

    @pMin.setter
    def pMin(self, newParam):
        #value in cliprange will not exceed available psi.
        self.pClipRange[0] = np.maximum(newParam, self.p.min())

    @property
    def pMax(self):
        return np.minimum(self.pClipRange[1], self.p.max())

    @pMax.setter
    def psiMax(self, newParam):
        #value in cliprange will not exceed available q.
        self.pClipRange[1] = np.minimum(newParam, self.p.max())

    @property
    def pLimsString(self):
        return u"{0:.3g} ≤ psi ({psiMagnitudeName}) ≤ {1:.3g}".format(
                self.pUnit.toDisplay(self.pMin),
                self.pUnit.toDisplay(self.pMax),
                pMagnitudeName = self.pUnit.displayMagnitudeName)

    @property
    def dataContent(self):
        """shows the content of the loaded data: Q, I, IErr, etc"""
        content = []
        if self.q is not None:
            content.append('Q')
        if self.i is not None:
            content.append('I')
        if self.hasError:
            content.append('IErr')
        if self.is2d:
            content.append('Psi')
        return ", ".join(content)

    @property
    def filename(self):
        return self._filename

    def setFilename(self, fn):
        """Stores the absolute path to this data file.
        Should be reviewed when data sets can be created from several files."""
        if fn is None or not os.path.isfile(fn):
            return
        self._filename = os.path.abspath(fn)

    @classmethod
    def load(cls, filename):
        """Factory method for creating SASData objects from file."""
        if not os.path.isfile(filename):
            logging.warning("File '{0}' does not exist!".format(filename))
            return
        
        logging.info("Loading '{0}' ...".format(filename))

        if str(filename[-4:]).lower() == '.pdh':
            sasFile = PDHFile(filename)
        else:
            sasFile = AsciiFile(filename) # works for CSV too

        sasData = cls(sasFile.name, sasFile.data)
        sasData.setFilename(sasFile.filename)
        return sasData

    def copy(self):
        cpy = SASData(self.title, self.origin)
        cpy.setFilename(self.filename)
        return cpy

    @property
    def q(self):
        """Q-Vector at which the intensities are measured."""
        return self.qUnit.toSi(self.origin[:, 0])

    @property
    def i(self):
        """Measured intensity at q."""
        return self.iUnit.toSi(self.origin[:, 1])

    @property
    def e(self):
        """Uncertainty or Error of the intensity at q loaded from file."""
        return self.iUnit.toSi(self.origin[:, 2])

    @property
    def u(self):
        """Corrected uncertainty or error of the intensity at q."""
        return self._uncertainty

    @property
    def p(self):
        """Psi-Vector."""
        return self.origin[:, 3] 

    @property
    def is2d(self):
        """Returns true if this dataset contains two-dimensional data with
        psi information available."""
        return self.origin.shape[1] > 3 # psi column is present

    @property
    def hasError(self):
        """Returns True if this data set has an error bar for its
        intensities."""
        return self.origin.shape[1] > 2

    @property
    def qUnit(self):
        return self._qUnit

    @property
    def pUnit(self):
        return self._pUnit

    @property
    def iUnit(self):
        return self._iUnit

    @property
    def rUnit(self):
        return self._rUnit

    def __init__(self, *args):
        #set unit definitions for display and internal units
        self._qUnit = ScatteringVector(u"nm⁻¹")
        self._pUnit = Angle(u"˚")
        self._iUnit = ScatteringIntensity(u"(m sr)⁻¹")
        self._rUnit = Length(u"nm")

        DataSet.__init__(self, *args)
        self._sizeEst = np.pi / np.array([self.q.max(),
                                          abs(self.q.min()) ])
        self._prepareUncertainty()

    def _prepareUncertainty(self):
        self._uncertainty = self.eMin * self.i
        minUncertaintyPercent = self.eMin * 100.
        if not self.hasError:
            logging.warning("No error column provided! Using {}% of intensity."
                            .format(minUncertaintyPercent))
        else:
            count = sum(self._uncertainty > self.e)
            if count > 0:
                logging.warning("Minimum uncertainty ({}% of intensity) set "
                                "for {} datapoints.".format(
                                minUncertaintyPercent, count))
            self._uncertainty = np.maximum(self._uncertainty, self.e)

    def __eq__(self, other):
        return (np.all(self.origin == other.origin)
                and self.title == other.title
                and self.filename == other.filename)

    def __neq__(self, other):
        return not self.__eq__(other)

    def sphericalSizeEst(self):
        return self._sizeEst

    def clip(self, *args, **kwargs):
        """Returns the clipped dataset."""
        clipped = self.copy()
        clipped.setOrigin(self.origin[self.clipMask(*args, **kwargs), :])
        return clipped

    def clipMask(self):
        #, qBounds = None, psiBounds = None,
        #               maskNegativeInt = False, maskZeroIntnt = False):
        """
        If q and/or psi limits are supplied in the dataset,
        clips the Dataset to within the set limits. 
        """
        # init indices: index array is more flexible than boolean masks
        validIndices = np.where(np.isfinite(self.q))[0]
        def cutIndices(validIndices, mask):
            validIndices = validIndices[mask]

        # Optional masking of negative intensity
        if self.maskZeroInt:
            # FIXME: compare with machine precision (EPS)?
            cutIndices(validIndices,self.i[validIndices] == 0.0)
        if self.maskNegativeInt:
            cutIndices(validIndices,self.i[validIndices] > 0.0)
        # clip to q bounds
        cutIndices(validIndices,self.q[validIndices] >= self.qMin)
        cutIndices(validIndices,self.q[validIndices] <= self.qMax)
        # clip to psi bounds
        if self.is2d:
            cutIndices(validIndices,self.p[validIndices] > self.pMin)
            cutIndices(validIndices,self.p[validIndices] <= self.pMax)

        return validIndices

if __name__ == "__main__":
    import doctest
    doctest.testmod()

# vim: set ts=4 sts=4 sw=4 tw=0:
