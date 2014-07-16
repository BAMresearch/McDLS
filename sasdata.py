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

import os # Miscellaneous operating system interfaces
import logging
import numpy as np # For arrays
from cutesnake.datafile import PDHFile, AsciiFile
from cutesnake.dataset import DataSet, DisplayMixin
from cutesnake.utils import isList
from cutesnake.utilsgui import processEventLoop

class SASData(DataSet, DisplayMixin):
    """Represents one set of data from a unique source (a file, for example).
    """
    _sizeEst = None
    _emin = 0.01 # minimum possible error (1%)
    _uncertainty = None
    _filename = None
    _qUnits = 1.e9 #can be set to scale to standard units (m^-1)
    _iUnits = 1. #can be set to scale to standard units ((m sr)^-1)
    _qClipRange = [-np.inf, np.inf] #manually set Q clip range for this dataset

    @staticmethod
    def displayDataDescr():
        return ("filename", "data points", "data content", 
                "Q limits", "est. sphere size")

    @property
    def displayData(self):
        return ("title", "count", "dataContent", 
                "qLimsString", "sphericalSizeEstText")

    @property
    def count(self):
        return len(self.q)

    @property
    def sphericalSizeEstText(self):
        return "min: {0:.4g}, max: {1:.4g}".format(*self.sphericalSizeEst())

    @property
    def qClipRange(self):
        return self._qClipRange

    def qMin(self):
        #returns minimum q from data or qClipRange, whichever is larger
        return np.maximum(self.qClipRange[0], self.q.min())

    def qMax(self):
        #returns maximum q from data or qClipRange, whichever is smaller
        return np.minimum(self.qClipRange[1], self.q.max())

    @property
    def qLimsString(self):
        return "min Q: {0:.4g}, max Q: {1:.4g}".format(
                self.qMin(), self.qMax())

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
    def qUnits(self):
        """Scaling factor for q to move to m^-1."""
        return self._qUnits

    @property
    def iUnits(self):
        """Scaling factor for q to move to m^-1."""
        return self._iUnits

    @property
    def q(self):
        """Q-Vector at which the intensities are measured."""
        return self.origin[:, 0] 

    @property
    def i(self):
        """Measured intensity at q."""
        return self.origin[:, 1]

    @property
    def e(self):
        """Uncertainty or Error of the intensity at q."""
        return self.origin[:, 2] * self.iUnits

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

    def __init__(self, *args):
        DataSet.__init__(self, *args)
        self._sizeEst = np.pi / np.array([self.q.max(),
                                          abs(self.q.min()) ])
        self._prepareUncertainty()

    def _prepareUncertainty(self):
        self._uncertainty = self.minUncertainty() * self.i
        minUncertaintyPercent = self.minUncertainty() * 100.
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

    def minUncertainty(self):
        return self._emin

    @property
    def uncertainty(self):
        return self._uncertainty

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
        clipped.setOrigin(self.origin[
                self.clipMask(*args, **kwargs)])
        return clipped

    def clipMask(self, qBounds = None, psiBounds = None,
                       maskNegativeInt = False, maskZeroInt = False):
        """If q and/or psi limits are supplied in self.Parameters,
        clips the Dataset to within the supplied limits. Copies data to
        :py:attr:`self.FitData` if no limits are set.
        """
        # init indices: index array is more flexible than boolean masks
        validIndices = np.where(np.isfinite(self.q))[0]
        def cutIndices(validIndices, mask):
            validIndices = validIndices[mask]

        # Optional masking of negative intensity
        if maskZeroInt:
            # FIXME: compare with machine precision (EPS)?
            cutIndices(validIndices,self.i[validIndices] == 0.0)
        if maskNegativeInt:
            cutIndices(validIndices,self.i[validIndices] > 0.0)
        if isList(qBounds):
            # excluding the lower q limit may prevent q = 0 from appearing
            cutIndices(validIndices,self.q[validIndices] > min(qBounds))
            cutIndices(validIndices,self.q[validIndices] <= max(qBounds))
        if isList(psiBounds) and self.is2d:
            cutIndices(validIndices,self.p[validIndices] > min(psiBounds))
            cutIndices(validIndices,self.p[validIndices] <= max(psiBounds))

        return validIndices

if __name__ == "__main__":
    import doctest
    doctest.testmod()

# vim: set ts=4 sts=4 sw=4 tw=0:
