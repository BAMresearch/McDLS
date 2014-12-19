# -*- coding: utf-8 -*-
# bases/datafile/asciifile.py

from numpy import array as np_array
from bases.datafile import DataFile
from utils.error import FileError
from utils import isString

class AsciiFile(DataFile):
    """A generic ascii data file."""
    valueFormat = "{0: 14.6E}" # format of data values for ascii export
    separator = " "
    newline="\n"

    @staticmethod
    def sanitizeData(data):
        data = DataFile.sanitizeData(data)
        assert data.ndim >= 2 # check for at least 2 dimensions
        return data

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
        with open(filename, mode) as fd:
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

    @classmethod
    def readRow(cls, fields, dataType = float, **kwargs):
        """Converts each field to the requested datatype.
           Raises an error if it is incompatible,
           the line is skipped in that case."""
        try:
            # just converted to tuple
            return tuple((dataType(f) for f in fields))
        except:
            raise ValueError

    @classmethod
    def readFile(cls, filename, dataType = float):
        fileData = []
        with open(filename, 'rU') as fd:
            linenr = 0
            for line in fd:
                linenr += 1
                # strip trailing white space, replace decimal operators 
                # eventually, split data fields
                # we read floating point numbers only
                fields = (line.strip()
                              .replace(",",".")
                              .replace(";"," ")
                              .split())
                record = None
                try:
                    # may raise exception
                    record = cls.readRow(fields, filename = filename,
                                         lineNumber = linenr,
                                         dataType = dataType)
                except ValueError:
                    continue
                if record is None or len(record) <= 0:
                    continue
                fileData.append(record)

        dataCount = len(fileData)
        if dataCount <= 0:
            raise FileError("No data columns found!", filename)
        return np_array(fileData, dataType)

# vim: set ts=4 sts=4 sw=4 tw=0: 
