# -*- coding: utf-8 -*-
# utils/hdf.py

"""
Helper functions for HDF5 functionality
"""

import logging
import inspect
import os.path
import h5py
from utils import isCallable, isString, isList, isNumber, isInteger, classname

# from utils.devtools import DBG

def getCallerInfo(referenceType = None, stackOffset = 0):
    """*referenceType*: Stop the search for a frame when this type for a local
    'self' is found.
    *stackOffset*: grab that frame counted from the last instead of search"""
    out = ""
    stack = inspect.stack()
    if isinstance(referenceType, type) and referenceType is not type(None):
        # search for the appropriate frame in user code
        frame = stack[stackOffset][0]
        while ('self' not in frame.f_locals
               or isinstance(frame.f_locals['self'], referenceType)):
            stackOffset += 1
            frame = stack[stackOffset][0]
    else:
        stackOffset += 2
    if len(stack) > stackOffset:
        frame = stack[stackOffset][0]
        head, fn = os.path.split(frame.f_code.co_filename)
#        func = inspect.getframeinfo(frame).function
        out = (u"({} l.{})".format(fn, frame.f_lineno))
    return out

class HDFWriter(object):
    """Represents an open HDF file location in memory and keeps track of
    the current address/name for reading or writing. Once this object looses
    scope, its data is actually written to file."""
    _handle = None
    _location = None

    def __init__(self, hdfHandle):
        self._handle = hdfHandle
        self._location = hdfHandle.name

    def __enter__(self):
        """Implements a *with* statement context manager."""
        self._handle.__enter__()
        return self

    def __exit__(self, *args, **kwargs):
        """Implements a *with* statement context manager."""
        return self._handle.__exit__(*args, **kwargs)

    @classmethod
    def open(cls, filename):
        return HDFWriter(h5py.File(filename, driver = 'core', backing_store = True))

    @property
    def location(self):
        return self._location

    def log(self, msg):
        logging.debug(u"[{}] {}".format(classname(self), msg))

    def _writeLocation(self):
        if self.location not in self._handle:
            self._handle.require_group(self.location)
        return self._handle[self.location]

    def writeAttributes(self, **kwargs):
        for key, value in kwargs.iteritems():
            # leading _ from keys are removed, they come from __getstate__()
            # which gathers member variables (usually private with leading _)
            self.writeAttribute(key.lstrip('_'), value)

    def writeAttribute(self, key, value):
        if value is None:
            return
        self.log("attribute '{loc}/{k}': '{v}'"
                 .format(loc = self.location.rstrip('/'), k = key, v = value))
        self._writeLocation().attrs[key] = value

    def writeDataset(self, name, data):
        if not isList(data):
            self.log("dataset {} is not of a list type! (={})"
                     .format(name, classname(data)))
            return
        shape = getattr(data, 'shape',
                        "(len: {})".format(len(data)))
        self.log("dataset '{loc}/{name}' {shape}"
                 .format(loc = self.location.rstrip('/'),
                         name = name, shape = shape))
        writeLocation = self._writeLocation()
#        DBG("loc: ", writeLocation)
        if name in writeLocation:
            del writeLocation[name]
        try:
            writeLocation.create_dataset(name, data = data,
                                         compression = "gzip")
        except TypeError:
            # a TypeError is raised for non-chunkable data (such as string)
            writeLocation.create_dataset(name, data = data)

    def writeMembers(self, obj, *members):
        assert(len(members))
        for member in members:
            self.writeMember(obj, member)

    def writeMember(self, obj, memberName):
        if isString(obj):
            logging.warning(u"String as object provided! "
                            + self._warningPrefix(obj, memberName))
        if isInteger(memberName) and isList(obj):
            member = obj[memberName]
            memberName = str(memberName)
        else:
            member = getattr(obj, memberName, None)
#        DBG(member)
        if member is None:
            self.log(u"skipped " + self._warningPrefix(obj, memberName)
                     + u"It is empty or does not exist (=None).")
            return
        if isCallable(member) and not isinstance(member, HDFMixin):
            # ParameterBase instances are callable but also HDFMixins
            member = member()
        if hasattr(member, "hdfWrite"): # support instances and types
            # store the member in a group of its own
            oldLocation = self.location
            self._location = "/".join((oldLocation.rstrip('/'), memberName))
            member.hdfWrite(self) # recursion entry, mind the loops!
            self._location = oldLocation
        elif isList(member):
            self.writeDataset(memberName, member)
        elif isString(member) or isNumber(member):
            self.writeAttribute(memberName, member)
        else:
            self.log(u"skipped " + self._warningPrefix(obj, memberName)
                     + "(={}) It is not a compatible value type!"
                        .format(classname(member)))

    def _warningPrefix(self, obj, memberName):
        return (u"{cls}.{mem} {cal} ".format(cal = getCallerInfo(type(self)),
                cls = classname(obj), mem = memberName))

class HDFMixin(object):

    def hdfStore(self, filename):
        """Writes itself to an HDF file at the given position or group."""
        with HDFWriter.open(filename) as hdf:
#            DBG(hdf)
            self._hdfWrite(hdf)

    def _hdfWrite(self, hdf):
        """A private wrapper for custom hdfWrite() methods in sub classes.
        Allows to handle common preprocessing without requiring to call the
        method in the parent class."""
        assert isinstance(hdf, HDFWriter)
#        DBG(hdf.location)
        self.hdfWrite(hdf)

    def hdfWrite(self, hdf):
        """To be overridden by sub classes to store themselves in an HDF
        structure. *hdf*: a HDFWriter instance."""
        logging.warning(u"Not defined: "
                        + hdf._warningPrefix(self, "hdfWrite()"))
        pass

    @classmethod
    def hdfLoad(self):
        """Restores an instance of this type from a given HDF file location
        or group."""
        pass

# vim: set ts=4 sts=4 sw=4 tw=0:
