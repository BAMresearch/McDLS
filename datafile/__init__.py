# -*- coding: utf-8 -*-
# bases/datafile/__init__.py

__all__ = ["DataFile", "AsciiFile", "ArrayFile", "PDHFile", "PDHHeader",
           "CGSFile", "getFileFilter", "loadfile"]

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
        for descr, ext in cls.extensions:
            filefilter.append(("{0} " + extFmt).format(descr, ext))
    import sys
    print >>sys.__stderr__, "filefilter", filefilter
    return filefilter

# vim: set ts=4 sw=4 sts=4 tw=0:
