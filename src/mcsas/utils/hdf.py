# -*- coding: utf-8 -*-
# utils/hdf.py

"""
Helper functions for HDF5 functionality
"""

import logging
import inspect
import os.path
import h5py
import numpy as np

from . import isCallable, isString, isList, isNumber, isInteger, classname

HDFString = str # Python 3 default
try:
    HDFString = unicode # fails with Python 3
except NameError:
    pass

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


def HDFCleanup(infile):
    """Unused space is reclaimed by making a copy of the contents in the 
    current hdf 5 file object, and moves the copy in place of the original. 
    if the input argument supplied is a file name instead of an HDF5 object, 
    the method returns nothing. Else, the method returns the new HDF5 object"""
    
    inputIsFilename = False
    if isString(infile): # a filename was supplied
        infile = h5py.File(infile)
        inputIsFilename = True
    def hdfDeepCopy(infile, outfilename):
        """Copies the internals of an open HDF5 file object (infile) to a second file."""
        infile.flush()
        outfile = h5py.File(outfilename, "w") # create file, truncate if exists
        for item in infile:
            outfile.copy(infile[item], item)
        for attr_name in infile.attrs:
            outfile.attrs[attr_name] = infile.attrs[attr_name]
        outfile.close()
    # shutil.move is more flexible than os.rename, and can move across disks.
    origfname = infile.filename
    # generate a temporary file for the original backup
    tempfname = "{}.hdf5".format(str(uuid.uuid4().hex))
    # temporary output filename
    tempofname = "{}.hdf5".format(str(uuid.uuid4().hex))
    shutil.copy(origfname, tempfname) # backup copy
    try: 
        hdfDeepCopy(infile, tempofname) # copy internals
        infile.close() # close old file
        shutil.move(tempofname, origfname) # replace with new
    except:
        # move back the backup
        shutil.move(tempfname, origfname) 
    infile = h5py.File(origfname) # reopen new file
    # cleanup:
    for filename in [tempfname, tempofname]:
        try:
            os.remove(filename)
        except OSError:
            pass
    if not inputIsFilename:
        return infile

class HDFWriter(object):
    """Represents an open HDF file location in memory and keeps track of
    the current address/name for reading or writing. Once this object looses
    scope, its data is actually written to file."""
    _handle = None
    _location = None

    def __init__(self, hdfHandle, rootLocation = None):
        self._handle = hdfHandle
        self._location = hdfHandle.name
        if rootLocation is not None:
            self._location = rootLocation

    def __enter__(self):
        """Implements a *with* statement context manager."""
        self._handle.__enter__()
        return self

    def __exit__(self, *args, **kwargs):
        """Implements a *with* statement context manager."""
        return self._handle.__exit__(*args, **kwargs)

    @classmethod
    def open(cls, filename, rootLocation = None):
        return HDFWriter(h5py.File(filename, 'w'), rootLocation)

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
        for key, value in kwargs.items():
            # leading _ from keys are removed, they come from __getstate__()
            # which gathers member variables (usually private with leading _)
            self.writeAttribute(key.lstrip('_'), value)

    def writeAttribute(self, key, value):
        if value is None:
            return
        if isString(value):
            value = HDFString(value)
        elif type(value) == type(True): # filter boolean
            value = np.int8(value)
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
        if member is None:
            self.log(u"skipped " + self._warningPrefix(obj, memberName)
                     + u"It is empty or does not exist (=None).")
            return

        if isCallable(member) and not hasattr(member, "hdfWrite"):
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

    def hdfStore(self, filename, rootLocation = None):
        """Writes itself to an HDF file at the given position or group."""
        with HDFWriter.open(filename, rootLocation) as hdf:
            self._hdfWrite(hdf)

    def _hdfWrite(self, hdf):
        """A private wrapper for custom hdfWrite() methods in sub classes.
        Allows to handle common preprocessing without requiring to call the
        method in the parent class."""
        assert isinstance(hdf, HDFWriter)
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
