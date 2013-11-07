# -*- coding: utf-8 -*-
# __init__.py

import sys
import collections
import platform
import hashlib
import os.path
import logging
import numpy

try:
    from cutesnake.qt import QtCore
    from QtCore import QString
except:
    # defines QString as regular python string if it could not be imported
    QString = basestring

EPS = sys.float_info.epsilon

# object tests

def isList(obj):
    return (not isString(obj) and
            (isinstance(obj, collections.Sequence)
             or (isinstance(obj, numpy.ndarray)
                 and obj.ndim < 2)))

def isString(obj):
    return isinstance(obj, basestring) or isinstance(obj, QString)

def isNonEmptyString(obj):
    return (isString(obj) and len(obj) > 0)

def isMap(obj):
    return isinstance(obj, collections.Mapping)

def isSet(obj):
    return isinstance(obj, collections.Set)

def isNumber(obj):
    # current implementation does not account for all possibilities.
    # trying to use numbers.Number and isinstance does not accept things
    # like float32:
    # return (isinstance(obj, int) or isinstance(obj, float) or
    #         isinstance(obj, long) or isinstance(obj, complex))

    # float(obj) gives false positive for strings like "1"
    # which are not supposed to be numbers
    if isString(obj):
        return False
    # the faster way to do this is "Duck typing", according to SO:
    try:
        float(obj)
    except StandardError: # catch ValueError and TypeError (for None)
        try:
            complex(obj)
        except StandardError:
            return False
    return True

def isInteger(obj):
    return (isinstance(obj, int) or isinstance(obj, long))

# environment tests

def isLinux():
    return platform.system().lower() in "linux"

def isMac():
    return platform.system().lower() in "darwin"

def isWindows():
    return platform.system().lower() in "windows"

def isFrozen():
    return hasattr(sys, "frozen")

# utility functions

def makeHash(s):
    m = hashlib.sha1()
    m.update(s)
    return m.hexdigest()

def makeFilename(path, basename, suffix, ext):
    """Creates a file name based on provided information."""
    return os.path.join(path, basename +'_'+ suffix +'.'+ ext)

def formatList(lst, precision = 4):
    fmt = ":." + str(precision) + "f"
    fmt = "; ".join(["{" + str(i) + fmt + "}"
        for i in range(0, len(lst))])
    return fmt.format(*lst)

#import inspect

def getTypeFromModules(name, modNames):
#    print "searching for", name
#    for i in range(0, 4):
#        print i
#        frame = inspect.stack()[i]
#        print "TEST", frame[0].f_globals.get('__name__', None)
    if not isList(modNames):
        modNames = (modNames, )
    objType = None
    for modName in modNames:
        try:
            mod = __import__(modName, globals(), locals(), (name, ))
        except StandardError:
            continue
        objType = getattr(mod, name, None)
        if objType is not None:
            break
    assert objType is not None, ("Object type {0} not found in modules {1}!"
                                 .format(name, modNames))
    return objType

def transformToString(t):
    """
    Generates a printable string in matrix form of a QTransform.

    >>> from utilsgui import printQTransform
    >>> printQTransform(None)
    ''
    >>> from PyQt4.QtGui import QTransform
    >>> printQTransform(  QTransform.fromScale(1.2, 0.7)
    ...                 * QTransform.fromTranslate(0.3, 1.4))
    '\\n+1.200e+00 +0.000e+00 +0.000e+00 \\n+0.000e+00 +7.000e-01 +0.000e+00 \\n+3.000e-01 +1.400e+00 +1.000e+00 '
    """
    res = ""
    if t is None:
        return res
    for row in [[t.m11(), t.m12(), t.m13()], \
                [t.m21(), t.m22(), t.m23()], \
                [t.m31(), t.m32(), t.m33()]]:
        res += "\n"
        for col in row:
            res += "{0:>+10.3e} ".format(col)
    return res

# some utility functions for testing

def log2file(text):
    """Logs some text to a temporary file.
    For debugging logging and output redirection and formatting."""
    with open("/tmp/test", "a") as fd:
        fd.write(bytearray(text, "utf8")+"\n")

def writeout(fn, arr):
    numpy.set_printoptions(threshold = 1e6)
    if isinstance(arr, (tuple, list)):
        arr = numpy.hstack([a.flatten() for a in arr])
    open(fn, 'w').write(repr(arr.flatten()))
    print "writing done"

def writecsv(fn, arr):
    # gnuplot -p -e 'set datafile separator ";"; set logscale xy; ...
    #     plot "/tmp/test.csv" using 1:2 with lines'
    import csv
    writer = csv.writer(open(fn, 'w'), delimiter=';')
    if len(arr.shape) == 1:
        arr = arr.reshape((-1, 1))
    for line in tuple(arr):
        writer.writerow(line)
    print "writing done"
    #No closing of the file?

# vim: set ts=4 sw=4 sts=4 tw=0:
