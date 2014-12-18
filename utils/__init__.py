# -*- coding: utf-8 -*-
# __init__.py

import numpy as np
from cutesnake.utils import (isList, isString, isNonEmptyString, isMap, isSet,
                             isNumber, isInteger)
from cutesnake.utils import isLinux, isMac, isWindows, isFrozen
from cutesnake.utils import testfor, assertName
from cutesnake.utils import mixedmethod
from cutesnake.utilsgui import processEventLoop

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

# vim: set ts=4 sw=4 sts=4 tw=0:
