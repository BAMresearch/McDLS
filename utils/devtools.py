# -*- coding: utf-8 -*-
# utils/devtools.py

"""
Definitions used during development only.
This module is supposed to not being required in a release package.
"""

from __future__ import print_function
from builtins import str
import sys
import inspect
import os, os.path
import tempfile
import codecs

prefix = os.path.splitext(os.path.basename(__file__))[0] + "_"
_logfd, _logfn = tempfile.mkstemp(prefix = prefix)
os.close(_logfd) # not used/needed
print("DBGF() to: '{}'".format(_logfn), file = sys.__stderr__)

def DBGF(*args):
    with codecs.open(_logfn, 'a', encoding = 'utf8') as fd:
        fd.write(_formatMessage(*args) + u"\n")

def _messagePrefix():
    stack = inspect.stack()
    res = []
    if len(stack) > 1:
        frame = stack[3][0]
        head, fn = os.path.split(frame.f_code.co_filename)
        head, mod = os.path.split(head)
        fn = os.path.splitext(fn)[0]
        func = inspect.getframeinfo(frame).function
        prefix = u"[{0:04d}|{1}]".format(frame.f_lineno, '.'.join((mod, fn)))
        res += [func, prefix]
    return res

def _formatMessage(*args):
    return u" ".join(_messagePrefix() + [str(a) for a in args])

def DBG(*args, **kwargs):
    _file = kwargs.pop('_file', sys.__stderr__)
    print(_formatMessage(*args), file = _file)

# vim: set ts=4 sts=4 sw=4 tw=0:
