# -*- coding: utf-8 -*-
# datafile/datafile.py

from builtins import str
from builtins import object
import os.path
from abc import ABCMeta, abstractmethod, abstractproperty
from numpy import ndarray as np_ndarray
from future.utils import with_metaclass

from ..utils.error import FileError
from ..utils.lastpath import LastPath
from ..utils import classproperty

class DataFile(with_metaclass(ABCMeta, object)):
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
    _filename = None

    @abstractproperty
    @classproperty
    @classmethod
    def fileFilter(cls):
        """Extensions of this file type. No extension by default.
           Use it to construct a proper file name. Reimplement in subclasses.
           """
        raise NotImplementedError

    @classproperty
    @classmethod
    def extensions(cls):
        return (e for t, e in cls.fileFilter)

    @property
    def filename(self):
        """Absolute path name of this file."""
        return self._filename

    def setFilename(self, filename):
        """Checks provided filename for plausibility and updates LastPath"""
        self._filename = self.sanitizeReadFilename(filename)
        LastPath.set(self.filename)

    @staticmethod
    def sanitizeReadFilename(filename):
        """Checks provided filename for plausibility and updates LastPath."""
        filename = os.path.abspath(str(filename))
        if not os.path.isfile(filename):
            raise FileError("Given file does not exist!", filename)
        return filename

    @property
    def name(self):
        """The plain name of the file with path and extension stripped."""
        return os.path.basename(os.path.splitext(self.filename)[0])

    # helpers for reading

    @abstractmethod
    def readFile(self, **kwargs):
        """Gets a proper file name and returns file data.
        May modify the instance. To be reimplemented."""
        raise NotImplementedError

    @abstractmethod
    def getDataObj(self):
        """Creates and returns the appropriate DataObj instance for this file
        type."""
        raise NotImplementedError

    def __init__(self, filename, **kwargs):
        self.setFilename(filename)
        self.readFile(**kwargs)

    # helpers for writing

    @classmethod
    def sanitizeWriteFilename(cls, filename):
        """Checks and sets the file name to write to."""
        if len(cls.extensions) > 0:
            dummy, ex = cls.extensions[0]
            if filename[-len(ex):] != ex:
                filename = "{0}.{1}".format(filename, ex)
        filename = os.path.abspath(str(filename))
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
        self.setFilename(filename)

    @classmethod
    @abstractmethod
    def writeFile(cls, filename, data, **kwargs):
        """Gets a proper file name and numpy array and writes it to file.
        Reimplement this."""
        raise NotImplementedError

# vim: set ts=4 sts=4 sw=4 tw=0: 
