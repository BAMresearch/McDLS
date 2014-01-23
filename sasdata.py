# -*- coding: utf-8 -*-
# mcsas/sasdata.py

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
    _prepared = None

    @staticmethod
    def displayDataDescr():
        return ("filename", "data points", "est. sphere size")

    @property
    def displayData(self):
        return ("title", "count", "sphericalSizeEstText")

    @property
    def count(self):
        return len(self.q())

    @property
    def sphericalSizeEstText(self):
        return "min: {0:.4f}, max: {1:.4f}".format(*self.sphericalSizeEst())

    @property
    def filename(self):
        return self._filename

    def setFilename(self, fn):
        """Stores the absolute path to this data file.
        Should be reviewed when data sets can be created from several files."""
        if not os.path.isfile(fn):
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

    def q(self):
        return self.origin[:, 0]

    def i(self):
        return self.origin[:, 1]

    @property
    def is2d(self):
        """Returns true if this dataset contains two-dimensional data with
        psi information available."""
        return self.origin.shape[1] > 3 # psi column is present

    @property
    def prepared(self):
        return self._prepared

    @prepared.setter
    def prepared(self, data):
        self._prepared = data

    def __init__(self, *args):
        DataSet.__init__(self, *args)
        self._sizeEst = np.pi / np.array([self.q().max(),
                                          min(abs(self.q().min()),
                                              abs(np.diff(self.q()).min()))
                                         ])
        self.prepareUncertainty()

    def sphericalSizeEst(self):
        return self._sizeEst

    def minUncertainty(self):
        return self._emin

    def prepareUncertainty(self):
        self._uncertainty = self.minUncertainty() * self.i()
        minUncertaintyPercent = self.minUncertainty() * 100.
        if self.origin.shape[1] < 3:
            logging.warning("No error column provided! Using {}% of intensity."
                            .format(minUncertaintyPercent))
        else:
            count = sum(self._uncertainty > self.origin[:, 2])
            if count > 0:
                logging.warning("Minimum uncertainty ({}% of intensity) set "
                                "for {} datapoints.".format(
                                minUncertaintyPercent, count))
            self._uncertainty = np.maximum(self._uncertainty, self.origin[:, 2])

    def uncertainty(self):
        return self._uncertainty

    def clip(self, *args, **kwargs):
        self.prepared = self.origin[
                self.clipMask(self.origin, *args, **kwargs)]

    @staticmethod
    def clipMask(data, qBounds = None, psiBounds = None,
                        maskNegativeInt = False, maskZeroInt = False):
        """If q and/or psi limits are supplied in self.Parameters,
        clips the Dataset to within the supplied limits. Copies data to
        :py:attr:`self.FitData` if no limits are set.
        """
        # some shortcut functions, not performance critical as this function
        # is called rarely
        def q(indices):
            return data[indices, 0]
        def intensity(indices):
            return data[indices, 1]
        def psi(indices):
            return data[indices, 3]

        # init indices: index array is more flexible than boolean masks
        validIndices = np.where(np.isfinite(data[:, 0]))[0]
        def cutIndices(mask):
            validIndices = validIndices[mask]

        # Optional masking of negative intensity
        if maskZeroInt:
            # FIXME: compare with machine precision (EPS)?
            cutIndices(intensity(validIndices) == 0.0)
        if maskNegativeInt:
            cutIndices(intensity(validIndices) > 0.0)
        if isList(qBounds):
            # excluding the lower q limit may prevent q = 0 from appearing
            cutIndices(q(validIndices) > min(qBounds))
            cutIndices(q(validIndices) <= max(qBounds))
        if isList(psiBounds) and data.shape[1] > 3: # psi in last column
            # excluding the lower q limit may prevent q = 0 from appearing
            cutIndices(psi(validIndices) > min(psiBounds))
            cutIndices(psi(validIndices) <= max(psiBounds))

        return validIndices

# vim: set ts=4 sts=4 sw=4 tw=0:
