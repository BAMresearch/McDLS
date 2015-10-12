# -*- coding: utf-8 -*-
# bases/datafile/__init__.py

# exports
from datafile import DataFile
from asciifile import AsciiFile
from pdhfile import PDHFile, PDHHeader
from cgsfile import CGSFile

# implementation
from utils import isLinux

def getFileFilter():
    extFmt = "(*.{1})"
    if isLinux():
        # the extension in parentheses is not shown on linux, add it here
        extFmt = "*.{1} (*.{1})"
    filefilter = []
    for cls in AsciiFile, PDHFile, CGSFile:
        for descr, ext in cls.extensions:
            filefilter.append(("{0} " + extFmt).format(descr, ext))
    import sys
    print >>sys.__stderr__, "filefilter", filefilter
    return filefilter

# vim: set ts=4 sw=4 sts=4 tw=0:
