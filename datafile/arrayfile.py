# -*- coding: utf-8 -*-
# datafile/arrayfile.py

from __future__ import absolute_import # PEP328
from numpy import array as np_array
from numpy import ndarray as np_ndarray
from datafile.asciifile import AsciiFile
from utils.classproperty import classproperty
from dataobj.sasdata import SASData

class ArrayFile(AsciiFile):
    """A data file containing a single array of data, mostly."""
    _rawArray = None

    @classproperty
    @classmethod
    def fileFilter(cls):
        return (("All files", "*"),
                ("ASCII single-array", "dat"),)

    # helpers for reading

    @property
    def rawArray(self):
        return self._rawArray

    @rawArray.setter
    def rawArray(self, rawArray):
        """Sets a numpy.array of raw file data as instance attribute."""
        assert isinstance(rawArray, np_ndarray)
        assert rawArray.ndim >= 2 # check for at least 2 dimensions
        self._rawArray = rawArray

    def parseLines(self, asciiLines, **kwargs):
        """Parses lines of an ASCII file in order to extract a single array
        of numbers. Reimplement this in subclasses for different behaviour.
        """
        lastLine, rawArray = self.readArray(asciiLines, **kwargs)
        self.rawArray = rawArray

    def getDataObj(self):
        sasData = SASData(title = self.name, rawArray = self.rawArray)
        sasData.setFilename(self.filename)
        return sasData

# vim: set ts=4 sts=4 sw=4 tw=0: 
