# -*- coding: utf-8 -*-
# utils/devtools.py

"""
Definitions used during development only.
This module is supposed to not being required in a release package.
"""

from __future__ import print_function
import sys
import inspect
import os.path

def DBG(*args):
    stack = inspect.stack()
    prefix = ""
    if len(stack) > 1:
        frame = stack[1][0]
        head, fn = os.path.split(frame.f_code.co_filename)
        head, mod = os.path.split(head)
        fn = os.path.splitext(fn)[0]
        func = inspect.getframeinfo(frame).function
        prefix = u"[{0:04d}|{1}]".format(frame.f_lineno, '.'.join((mod, fn)))
    print(u" ".join([prefix, func]
        + [unicode(a) for a in args]), file = sys.__stderr__)

# vim: set ts=4 sts=4 sw=4 tw=0:
