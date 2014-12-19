# -*- coding: utf-8 -*-
# bases/datafile/datafile.py

import os.path
from abc import ABCMeta, abstractmethod
from numpy import ndarray as np_ndarray
from utils.error import FileError
from utils.lastpath import LastPath

class DataFile(object):
    """
    Base class for handling data files.
    Can be initialized with a file (name) to read or with a data array.

    Test error behaviour
    >>> from utils import DataFile, getTempFileName, getTempFile, FileNotFound
    >>> fn = getTempFileName()
    >>> try: DataFile.checkFilename(fn)
    ... except FileNotFound, e: str(e).find(fn) > 0
    True

    Prepare test data file
    >>> fd = getTempFile()
    >>> l = ['123 234\\n', '1,23 43.4\\r\\n', '2.3; 34,4\\n', 
    ...      '21.2 42 2\\n', '23,2 3.4 \\n']
    >>> fd.writelines(l)
    >>> fd.close()

    Test data parsing
    >>> df = DataFile()
    >>> df.loadFile(fd.name)
    [('123', '234'), ('1.23', '43.4'), ('2.3', '34.4'), ('21.2', '42', '2'), ('23.2', '3.4')]

    Remove test data finally
    >>> import os
    >>> if os.path.isfile(fd.name):
    ...     os.remove(fd.name)
    """
    __metaclass__ = ABCMeta
    _filename = None
    _data = None
    """Extensions of this file type. Starts with default extension.
       Use it to construct a proper file name."""
    extensions = ()

    @property
    def filename(self):
        """Absolute path name of this file."""
        return self._filename

    @filename.setter
    def filename(self, filename):
        """Checks provided filename for plausibility and updates LastPath"""
        self._filename = self.sanitizeReadFilename(filename)
        LastPath.set(self.filename)

    @staticmethod
    def sanitizeReadFilename(filename):
        """Checks provided filename for plausibility and updates LastPath."""
        filename = os.path.abspath(unicode(filename))
        if not os.path.isfile(filename):
            raise FileError("Given file does not exist!", filename)
        return filename

    @property
    def name(self):
        """The plain name of the file with path and extension stripped."""
        return os.path.basename(os.path.splitext(self.filename)[0])

    @property
    def data(self):
        """Numpy array of raw file data."""
        return self._data

    @data.setter
    def data(self, data):
        """Expects a numpy array."""
        self._data = self.sanitizeData(data)

    @staticmethod
    def sanitizeData(data):
        assert isinstance(data, np_ndarray)
        return data

    @property
    def hasData(self):
        return isinstance(self._data, np_ndarray)

    def __init__(self, filenameOrArray, **kwargs):
        try:
            self.data = filenameOrArray
        except StandardError, e:
            self.filename = filenameOrArray
            self.data = self.readFile(self.filename, **kwargs)

    @classmethod
    @abstractmethod
    def readFile(cls, filename, **kwargs):
        """Gets a proper file name and returns a numpy array.
        Reimplement this."""
        raise NotImplementedError

    @classmethod
    def sanitizeWriteFilename(cls, filename):
        """Checks and sets the file name to write to."""
        if len(cls.extensions) > 0:
            ex = cls.extensions[0]
            if filename[-len(ex):] != ex:
                filename = "{0}.{1}".format(filename, ex)
        filename = os.path.abspath(unicode(filename))
        dirname = os.path.dirname(filename)
        if not os.path.isdir(dirname):
            raise FileError("Directory does not exist:", dirname)
        return filename

    @classmethod
    def writeData(cls, filename, data, **kwargs):
        """Convenience method to write data directly to file."""
        self = cls(data, **kwargs)
        self.write(filename, **kwargs)

    def write(self, filename, **kwargs):
        filename = self.sanitizeWriteFilename(filename)
        self.writeFile(filename, self.data, **kwargs)
        self.filename = filename

    @classmethod
    @abstractmethod
    def writeFile(cls, filename, data, **kwargs):
        """Gets a proper file name and numpy array and writes it to file.
        Reimplement this."""
        raise NotImplementedError

# vim: set ts=4 sts=4 sw=4 tw=0: 
