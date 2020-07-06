# -*- coding: utf-8 -*-
# utils/__init__.py

from builtins import str
import sys
import codecs
import numpy as np
from .tests import (isList, isString, isNonEmptyString, isMap, isSet,
                         isNumber, isInteger, isCallable)
from .tests import isLinux, isMac, isWindows, isFrozen
from .tests import testfor, assertName
from .mixedmethod import mixedmethod
from .classproperty import classproperty

EPS = sys.float_info.epsilon

def clip(value, minv, maxv):
    """Expects a range tuple or list consisting of lower and upper limits."""
    if minv is None or maxv is None:
        return value
    # clips value to within set min/max limits.
    valueType = type(value)
    # clip to min/max independent on if value is list, int, float or array
    # print "arrayval: {}, min: {}, max: {}".format(value, selforcls.min(), selforcls.max())
    value = np.clip(np.array(value), minv, maxv)
    # return to original type:
    value = valueType(value)
    return value

def hashNumpyArray(arr):
    # get a read-only view
    view = arr.view()
    view.setflags(write = False)
    try:
        return hash(view.data)
    except ValueError:
        # Python 3: memoryview: hashing is restricted to formats 'B', 'b' or 'c'
        return hash(view.data.tobytes()) # python3

# fixes IOError on Windows: Errno 2, no such fn or dir
# http://stackoverflow.com/a/11442604
# better use str.encode(sys.getfilesystemencoding())
# and   os.cwd().decode(sys.getfilesystemencoding())
# even better: In Windows the default total path length must not
# exceed 260 characters
# (drive + :\ + 256 characters of filename + null terminator)
# http://superuser.com/a/811155

_winFilenamePrefix = u"\\\\?\\"
def fixFilename(filename):
    """Works around Windows file path length limitation of 260 chars."""
    if not isWindows() or filename.startswith(_winFilenamePrefix):
        return filename
    # work around maxlen=260 chars, \\?\ allows 32k chars max
    return _winFilenamePrefix + str(filename)

def mcopen(fn, mode, encoding = "utf8"):
    if 'b' in mode:
        return open(fixFilename(fn), mode)
    else:
        return codecs.open(fixFilename(fn), mode, encoding = encoding)

def classname(obj):
    if not isinstance(obj, type):
        obj = type(obj)
    return obj.__name__

# vim: set ts=4 sw=4 sts=4 tw=0:
