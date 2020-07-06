# -*- coding: utf-8 -*-
# datafile/asciifile.py

import codecs
from abc import ABCMeta, abstractmethod
from numpy import array as np_array
from datafile import DataFile
from utils.error import FileError
from utils import isString, isWindows, mcopen
from future.utils import with_metaclass

class AsciiFile(with_metaclass(ABCMeta, DataFile)):
    """A generic ascii data file."""
    valueFormat = "{0: 14.6E}" # format of data values for ascii export
    separator = " "
    newline="\n"

    # helpers for writing

    @classmethod
    def formatValue(cls, value):
        try:
            return cls.valueFormat.format(value)
        except ValueError:
            return "{0}".format(value)

    @classmethod
    def formatRow(cls, row, **kwargs):
        return cls.separator.join([cls.formatValue(value) for value in row])

    @classmethod
    def formatData(cls, data, **kwargs):
        return cls.newline.join([cls.formatRow(row, **kwargs)
                                 for row in data])

    @staticmethod
    def _write(filename, mode, asciiData):
        with mcopen(filename, mode) as fd:
            fd.write(asciiData)

    @classmethod
    def writeFile(cls, filename, data, **kwargs):
        cls._write(filename, 'w', cls.formatData(data, **kwargs))

    @classmethod
    def appendFile(cls, filename, data, **kwargs):
        """like writeFile but appends data to an existing file"""
        cls._write(filename, 'a', cls.formatData(data, **kwargs))

    @classmethod
    def _formatHeader(cls, header):
        if not isString(header):
            header = cls.separator.join(header)
        header += cls.newline
        return header

    @classmethod
    def writeHeaderLine(cls, filename, header):
        """writes a single-line header to a file consisting of a string or
        tuple of strings to be joined"""
        cls._write(filename, 'w', cls._formatHeader(header))

    @classmethod
    def appendHeaderLine(cls, filename, header):
        """writes a single-line header to a file consisting of a string or
        tuple of strings to be joined"""
        cls._write(filename, 'a', cls._formatHeader(header))

    # helpers for reading

    def readTuple(self, fields, dataType = float, **kwargs):
        """Converts each field to the requested datatype.
           Raises an error if it is incompatible,
           the line is skipped in that case."""
        try:
            # just converted to tuple
            return tuple((dataType(f) for f in fields))
        except:
            raise ValueError

    def readFile(self, **kwargs):
        asciiLines = None
        try:
            with mcopen(self.filename, 'r') as fd:
                asciiLines = fd.readlines()
        except UnicodeDecodeError:
            with mcopen(self.filename, 'r', encoding = 'latin1') as fd:
                asciiLines = fd.readlines()
        self.parseLines(asciiLines, **kwargs)

    @abstractmethod
    def parseLines(self, asciiLines, **kwargs):
        """Parses lines of an ASCII file in order to extract a single array
        of numbers. Reimplement this in subclasses for different behaviour.
        """
        raise NotImplementedError

    def readArray(self, asciiLines, dataType = float,
                  startLine = 0, endLine = None, **kwargs):
        """Reads a numpy.array from a specified segment (startLine, endLine)
        of a line buffer given by asciiLines. Stops at lines incompatible
        to previous lines read due to different number of fields or
        incompatible data type. Returns the last line successfully parsed and
        the populated numpy.array.
        """
        recordList = []
        for linenr, line in enumerate(asciiLines[startLine:endLine]):
            linenr += startLine
            # strip trailing white space, replace decimal operators 
            # eventually, split data fields
            # we read floating point numbers only
            if '.' in line:
                line = line.replace(","," ") # comma separated
            else: # no points in line
                # convert D/A/CH decimal separator to point
                line = line.replace(",",".")
            fields = (line.strip()
                          .replace(";"," ")
                          .split())
            record = None
            try:
                # may raise exception
                record = self.readTuple(fields, lineNumber = linenr,
                                        dataType = dataType)
                if not len(record): # ignore empty tuples
                    record = None
            except ValueError:
                pass # ignore it for now, record == None
            if record is None: # on parse failure of current line
                if not len(recordList):
                    continue # if still no compatible data found
                else:
                    break # data listing ends here
            elif len(recordList) and len(recordList[-1]) != len(record):
                break # do not append records of different size
            recordList.append(record)

        endLine = linenr - 1 # could not read last line
        recordCount = len(recordList)
        if recordCount <= 0:
            raise FileError("No data columns found!", self.filename)
        return endLine, np_array(recordList, dataType)

# vim: set ts=4 sts=4 sw=4 tw=0: 
