# -*- coding: utf-8 -*-
# tests.py

"""
Utils for testing something.
"""

from cutesnake.utils import isString

def testfor(condition, exception, errorMessage):
    if __debug__ and not condition:
        raise exception(errorMessage)

def assertName(newName, errorType, noWhitespace = False):
    testfor(isString(newName) and len(newName) > 0,
            errorType, "A name is mandatory!")
    testfor(not noWhitespace or newName.find(" ") < 0,
            errorType, "A name must not contain white space!")

# vim: set ts=4 sts=4 sw=4 tw=0:
