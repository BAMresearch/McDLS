# -*- coding: utf-8 -*-
# datafile/pdhfile.py

from __future__ import absolute_import # PEP328
from utils import isList, isInteger, isString
from utils.classproperty import classproperty
from datafile import AsciiFile

class PDHFile(AsciiFile):

    @classproperty
    @classmethod
    def extensions(cls):
        return (("PDH Saxsquant", "pdh"),)

    @classmethod
    def formatData(cls, data, description = None):
        hdr = PDHHeader(data.shape[0], description)
        asciiData = AsciiFile.formatData(data)
        asciiData = "{0}\n{1}".format(hdr, asciiData)
        return asciiData

    def parseLines(self, asciiLines, **kwargs):
        lastLine, rawArray = self.readArray(asciiLines, startLine = PDHHeader.maxLines, **kwargs)
        self.setRawArray(rawArray)

class PDHHeader(object):
    _entries = { # entry names and their header position (line, column)
        "description": (0, ),
        "keywords": (1, ),
        "dataCount": (2, 0),
        "sampleDetectorDistanceMM": (3, 1),
        "normalizationFactor": (3, 3),
        "waveLength": (3, 4),
        "dummy": (4, 0),
    }
    _values = None # list of lists

    @classmethod
    def _setup(cls):
        """Creates setter methods for each name in _entries."""
        def setValueFactory(index):
            return lambda self, value: self._setValue(index, value)
        for key, index in cls._entries.iteritems():
            key = "set" + key[0].upper() + key[1:]
            setattr(cls, key, setValueFactory(index))

    def _setValue(self, index, value):
        """Common setter method for all _entries."""
        assert isList(index)
        if not isList(self._values):
            self._values = []
        # enforce value types instead of testing it
        row = index[0]
        while len(self._values) <= row:
            self._values.append([])
        if len(index) == 1:
            # set the full line to value
            self._values[row] = value
        else:
            col = colCount = index[1]
            if row == 2:
                colCount = 8
            elif row > 2:
                colCount = 5
            line = self._values[row]
            while len(line) < colCount:
                line.append(0)
            line[col] = value

    @classproperty
    @classmethod
    def maxLines(cls):
        return max([index[0] for index in cls._entries.itervalues()]) + 1

    def __init__(self, dataCount, description = None):
        object.__init__(self)
        if any([not hasattr(self, key) for key in self._entries.iterkeys()]):
            self._setup()
        # init header with default values
        if not isString(description):
            description = ""
        self.setDescription(description)
        self.setKeywords(["SAXS", "BOX"])
        assert isInteger(dataCount)
        self.setDataCount(dataCount)
        self.setWaveLength(0) # setting spare lines
        self.setDummy(0)

    def line(self, index):
        """Returns the specified line of the header as string."""
        # 0 <= index <= self.lineCount
        index = max(index, 0)
        index = min(index, self.maxLines - 1)
        value = self._values[index]
        if index == 0:
            return "{0}".format(value)
        elif index == 1:
            return " ".join(value)
        elif index == 2:
            return " ".join(["{0: 9d}".format(v) for v in value])
        elif index > 2:
            return " ".join(["{0: 14.6E}".format(v) for v in value])

    def __str__(self):
        return "\n".join([self.line(i) for i in range(0, self.maxLines)])

# vim: set ts=4 sw=4 sts=4 tw=0:
