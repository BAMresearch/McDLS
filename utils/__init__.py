# -*- coding: utf-8 -*-
# __init__.py

from cutesnake.utils import (isList, isString, isNonEmptyString, isMap, isSet,
                             isNumber, isInteger)
from cutesnake.utils import isLinux, isMac, isWindows, isFrozen
from cutesnake.utils import testfor, assertName
from cutesnake.utils import mixedmethod
from cutesnake.utilsgui import processEventLoop

def clamp(value, inRange):
    """Expects a valid input range."""
    minmax = min(inRange), max(inRange)
    if value < minmax[0]:
        return minmax[0]
    elif value > minmax[1]:
        return minmax[1]
    return value

# vim: set ts=4 sw=4 sts=4 tw=0:
