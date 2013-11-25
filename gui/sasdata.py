# -*- coding: utf-8 -*-
# sasdata.py

import sys
import logging
import time
import os.path
import ConfigParser
import numpy
import numpy as np
import pickle
from cutesnake.qt import QtCore
from QtCore import QUrl
from cutesnake.dataset import DataSet, DisplayMixin
from cutesnake.utils import isList
from cutesnake.utils.lastpath import LastPath
from cutesnake.datafile import PDHFile, AsciiFile
from cutesnake.utilsgui.displayexception import DisplayException
from cutesnake.log import timestamp, addHandler
import cutesnake.log as log
from McSAS import McSAS

class SASData(DataSet, DisplayMixin):
    _sizeEst = None
    _emin = 0.01 # minimum possible error (1%)
    _uncertainty = None
    _filename = None

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

    @classmethod
    def load(cls, filename):
        if not os.path.isfile(filename):
            logging.warning("File '{0}' does not exist!".format(filename))
            return

        #probably not the right way of saving the filename somewhere:
        cls._filename = filename

        logging.info("Loading '{0}' ...".format(filename))

        if str(filename[-4:]).lower() == '.pdh':
            sasFile = PDHFile(filename)
        else:
            sasFile = AsciiFile(filename) # works for CSV too

        sasData = cls(sasFile.name, sasFile.data)
        return sasData

    def q(self):
        return self.origin[:, 0]

    def i(self):
        return self.origin[:, 1]

    def __init__(self, *args):
        DataSet.__init__(self, *args)
        self._sizeEst = np.array([np.pi/np.max(self.q()), np.pi/np.min(self.q())])
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

# vim: set ts=4 sts=4 sw=4 tw=0:
