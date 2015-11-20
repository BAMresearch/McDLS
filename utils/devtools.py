# -*- coding: utf-8 -*-
# utils/devtools.py

"""
Definitions used during development only.
This module is supposed to not being required in a release package.
"""

import sys

def DBG(*args):
    print >>sys.__stderr__, " ".join([unicode(a) for a in args])

# vim: set ts=4 sts=4 sw=4 tw=0:
