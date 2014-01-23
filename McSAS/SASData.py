# -*- coding: utf-8 -*-
# SASData.py

import numpy # For arrays
from numpy import (array, pi, diff)
import os # Miscellaneous operating system interfaces
from cutesnake.dataset import DataSet
from cutesnake.utils import isList
from cutesnake.utilsgui import processEventLoop

class SASData(DataSet):
    """Represents one set of data from a unique source (a file, for example).
    """
    _prepared = None
    _sizeBounds = None

    @classmethod
    def load(cls, filename):
        """Factory method for creating SASData objects from file."""
        if not os.path.isfile(filename):
            logging.warning("File '{0}' does not exist!".format(filename))
            return
        logging.info("Loading '{0}' ...".format(filename))
        sasFile = PDHFile(filename)
        sasData = cls(sasFile.name, sasFile.data)
        return sasData

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

    #this cannot be a @classmethod, apparently, as it will be a class ABC, noto DataSet
    def __init__(self, title, data):
        DataSet.__init__(self, title, data)

    #ibid.
    def setOrigin(self, data):
        DataSet.setOrigin(self, data)
        # determining sizeBounds from q vector
        q = data[:, 0]
        self._sizeBounds = pi / array([q.max(),
                                       min(abs(q.min()),
                                           abs(diff(q).min()))])

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
        validIndices = numpy.where(numpy.isfinite(data[:, 0]))[0]
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

    @property
    def sizeBounds(self):
        return self._sizeBounds

# vim: set ts=4 sts=4 sw=4 tw=0:
