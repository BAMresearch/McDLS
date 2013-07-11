# -*- coding: utf-8 -*-
# tests_test.py

"""
Testing test utils.
"""

from nose.tools import raises, nottest
from cutesnake.utils.tests import testfor, assertName

nottest(testfor)

class TestException(StandardError):
    pass

@raises(TestException)
def testtestfor():
    testfor(False, TestException, "test error")

@raises(TestException)
def testAssertNameNone():
    assertName(None, TestException)

@raises(TestException)
def testAssertNameEmpty():
    assertName("", TestException)

def testAssertNameOk():
    assertName("blub bla", TestException)

@raises(TestException)
def testAssertNameWhitespace():
    assertName("blub bla", TestException, noWhitespace = True)

def testAssertNameNoWhitespace():
    assertName("blubbla", TestException, noWhitespace = True)

# vim: set ts=4 sts=4 sw=4 tw=0:
