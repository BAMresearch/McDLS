# -*- coding: utf-8 -*-
# datafile/__init__.py

__all__ = ["DataFile", "AsciiFile", "ArrayFile", "PDHFile", "PDHHeader",
           "CGSFile", "getFileFilter", "loaddatafile"]

import logging
import os.path

from datafile import DataFile
from arrayfile import ArrayFile
from asciifile import AsciiFile
from pdhfile import PDHFile, PDHHeader
from cgsfile import CGSFile

from utils import isLinux, isString

def getFileFilter():
    """Returns the file filter text of all available data file formats which
    can be used with a file selection dialog UI."""
    extFmt = "(*.{1})"
    if isLinux():
        # the extension in parentheses is not shown on linux, add it here
        extFmt = "*.{1} (*.{1})"
    filefilter = []
    for cls in ArrayFile, PDHFile, CGSFile:
        for descr, ext in cls.fileFilter:
            filefilter.append(("{0} " + extFmt).format(descr, ext))
    return filefilter

def loaddatafile(filename):
    """Creates a DataFile instance from the given file name by selecting
    the appropriate specialised file object based on the file name extension
    or parts of the file contents."""

    if not isString(filename) or not os.path.isfile(filename):
        logging.warning("File '{0}' does not exist!".format(filename))
        return
    logging.info("Loading '{0}' ...".format(filename))
    path, ext = os.path.splitext(filename)
    ext = ext[1:].lower()
    datafile = None
    if ext in PDHFile.extensions:
        datafile = PDHFile(filename)
    elif ext in CGSFile.extensions:
        datafile = CGSFile(filename)
    else:
        datafile = ArrayFile(filename) # works for CSV too
    return datafile

# vim: set ts=4 sw=4 sts=4 tw=0:
