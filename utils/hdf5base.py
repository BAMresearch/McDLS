# -*- coding: utf-8 -*-
# utils/hdf5base.py

"""
Helper functions for HDF5 functionality
"""

import logging
import inspect
import os.path
import h5py
import numpy
from abc import ABCMeta
from utils import isCallable, isString, isList, classname

def h5w(wloc, field, hDat, hType = "dataset"):
    """
        writes dataset *hDat* to HDF5 location *wloc*, deleting if exists
        htype can be "dataset" or "attribute"
    """
    # remove old field, only removes link, does not reclaim!
    # ideally, new h5py file should be generated on end:
    # http://stackoverflow.com/questions/11194927/deleting-information-from-an-hdf5-file
    if "dataset" in hType:
        if field in wloc:
            del wloc[field]
        try:
            wloc.create_dataset(field, data = hDat, compression = "gzip")
        except TypeError:
            # a TypeError is raised for non-chunkable data (such as string)
            wloc.create_dataset(field, data = hDat)
        except:
            raise
    else:
        wloc.attrs[field] = hDat

from utils.devtools import DBG

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
        logging.debug(u"  [{}] {}".format(classname(self), msg))

    def writeAttribute(self, key, value):
        if value is None:
            return
        self.log("attribute '{loc}/{k}': '{v}'"
                 .format(loc = self.location.rstrip('/'), k = key, v = value))
        self._handle[self.location].attrs[key] = value

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
        writeLocation = self._handle[self.location]
        DBG("loc: ", writeLocation)
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
        member = getattr(obj, memberName, None)
        DBG(member)
        if member is None:
            self.log(u"skipped " + self._warningPrefix(obj, memberName)
                     + u"It is empty or does not exist (=None).")
            return
        if isCallable(member):
            member = member()
        if isList(member):
            self.writeDataset(memberName, member)
        elif isinstance(member, HDF5Mixin):
            # store the member in a group of its own
            oldLocation = self.location
            DBG("loc0: ", self.location)
            self._location = "/".join((oldLocation.rstrip('/'), memberName))
            DBG("loc1: ", self.location)
            self._handle.require_group(self.location)
            member.hdfWrite(self) # recursion entry, mind the loops!
            self._location = oldLocation
            DBG("loc3: ", self.location)
        else:
            self.log(u"skipped " + self._warningPrefix(obj, memberName)
                     + "(={}) It is neither of a list type nor a {}!"
                        .format(classname(member), classname(HDF5Mixin)))
            return

    def _warningPrefix(self, obj, memberName):
        return (u"{cls}.{mem} {cal}:\n".format(cal = getCallerInfo(type(self)),
                cls = classname(obj), mem = memberName))

class HDF5Mixin(object):
    __metaclass__ = ABCMeta
    _h5Datasets = []
    _h5Attrs = []
    _h5Callers = []
    _h5LocAdd = ''

    # mixin object (?) to add a write function to other classes
#    def __init__(self, **kwargs):
#        super(HDF5Mixin, self).__init__(**kwargs)

    def hdfStore(self, filename):
        """Writes itself to an HDF file at the given position or group."""
        with HDFWriter.open(filename) as hdf:
            DBG(hdf)
            self._hdfWrite(hdf)

    def _hdfWrite(self, hdf):
        """A private wrapper for custom hdfWrite() methods in sub classes.
        Allows to handle common preprocessing without requiring to call the
        method in the parent class."""
        assert isinstance(hdf, HDFWriter)
        DBG(hdf.location)
        self.hdfWrite(hdf)

    def hdfWrite(self, hdf):
        """To be overridden by sub classes to store themselves in an HDF
        structure. *hdf*: a HDFWriter instance."""
        pass

    @classmethod
    def hdfLoad(self):
        """Restores an instance of this type from a given HDF file location
        or group."""
        pass

    def writeHDF(self, filename, loc, item = None):
        """
        Writes the vector to an HDF5 output file *filename*, at location *loc*.
        This location should be e.g. "/mcentry01/[ sas | dls ]data01/x0".
        """
        loc += self._h5LocAdd
        if item is not None:
            loc += "/" + item
        # process datasets and attribute calls
        logging.debug("Writing at loc: {}".format(loc))
        with h5py.File(filename) as h5f:
            wloc = h5f.require_group(loc) # unicode's no problem
            for field in self._h5Datasets:
                hDat = getattr(self, field, None)
                if hDat is not None:
                    logging.debug("*** dataset name: {}"
                            .format(field))
                    h5w(wloc, field, hDat, hType = "dataset")
            for att in self._h5Attrs:
                logging.debug("*** attribute field: {}, data: {}"
                        .format(att, self.name))
                h5w(h5f[loc], att, self.name, hType = "attribute")

        # process custom callers
        for call in self._h5Callers:
            caller = getattr(self, call, None)
            if caller is None:
                # skip it
                continue
            callerWFunc = getattr(caller, "writeHDF", None)
            if callerWFunc is not None:
                logging.debug("Calling HDF5 writer for {}, at loc: {}"
                        .format(call, loc))
                callerWFunc(filename, loc, item = call)
            else:
                logging.warning("item {} does not have writeHDF functionality"
                        .format(call))

# vim: set ts=4 sts=4 sw=4 tw=0:
