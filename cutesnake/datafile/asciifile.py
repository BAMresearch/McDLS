# -*- coding: utf-8 -*-
# asciifile.py

from numpy import array as np_array
from cutesnake.datafile import DataFile
from cutesnake.utils.error import FileError
from cutesnake.utils.translate import tr

class AsciiFile(DataFile):
    """A generic ascii data file."""
    valueFormat = "{0: 14.6E}" # format of data values for ascii export

    @staticmethod
    def sanitizeData(data):
        data = DataFile.sanitizeData(data)
        assert data.ndim >= 2 # check for at least 2 dimensions
        return data

    @classmethod
    def formatRow(cls, row, **kwargs):
        return " ".join([cls.valueFormat.format(value) for value in row])

    @classmethod
    def formatData(cls, data, **kwargs):
        return "\n".join([cls.formatRow(row, **kwargs) for row in data])

    @classmethod
    def writeFile(cls, filename, data, **kwargs):
        asciiData = cls.formatData(data, **kwargs)
        with open(filename, 'w') as fd:
            fd.write(asciiData)

    @classmethod
    def readRow(cls, fields, **kwargs):
        return tuple(fields) # just converted to tuple

    @classmethod
    def readFile(cls, filename, dataType = None):
        fileData = []
        with open(filename) as fd:
            linenr = 1
            for line in fd:
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
                                         lineNumber = linenr)
                except ValueError:
                    continue
                else:
                    if record is None or len(record) <= 0:
                        continue
                    fileData.append(record)
                linenr += 1

        dataCount = len(fileData)
        if dataCount <= 0:
            raise FileError(tr("No data columns found!"), filename)
        return np_array(fileData, dataType)

# vim: set ts=4 sts=4 sw=4 tw=0: 
