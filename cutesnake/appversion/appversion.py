# -*- coding: utf-8 -*-
# appversion.py

import inspect
import os.path
import re
from cutesnake.utils import isNonEmptyString

class AppVersion(object):
    """
    Stores version meta data.
    """
    _name = None
    _number = None
    _organizationName = None
    _organizationDomain = None
    _defaultSettings = None

    def __init__(self, programName = None,
                       versionNumber = None,
                       organizationName = None,
                       organizationDomain = None,
                       defaultSettings = None):
        assert isNonEmptyString(programName)
        self._name = programName
        assert isNonEmptyString(versionNumber)
        self._number = versionNumber
        if organizationName is not None:
            assert isNonEmptyString(organizationName)
            self._organizationName = organizationName
        if organizationDomain is not None:
            assert isNonEmptyString(organizationDomain)
            self._organizationDomain = organizationDomain
        if defaultSettings is not None:
            assert type(defaultSettings) is dict
            assert len(defaultSettings) > 0
            self._defaultSettings = defaultSettings

    def name(self):
        return self._name

    def number(self):
        # pywin32 has a problem with version nr like "20111110"
        return self._number

    def organizationName(self):
        return self._organizationName

    def organizationDomain(self):
        return self._organizationDomain

    def defaultSettings(self):
        return self._defaultSettings

    @classmethod
    def isValid(cls, other):
        return issubclass(type(other), cls)

    @staticmethod
    def updateFile(module, newversion):
        """Updates the version string within a given module.
        Replaces the version string in the source code textually.
        Assumes there is an AppVersion constructor call in the module.
        """
        fn = os.path.abspath(inspect.getfile(module))
        if not os.path.exists(fn):
            return
        if os.path.splitext(fn)[-1].endswith("pyc"):
            fn = fn[:-1]
        # read current file content
        content = None
        with open(fn) as fd:
            content = fd.read()
        # get version number, parse and update it
        match = re.search(r'(versionNumber\s+=\s+(\"[^\"]*\"))', content)
        oldversion = match.group(2)
        # replace version string and write back
        newversion = '"{0}"'.format(newversion)
        content = content.replace(oldversion, newversion)
        with open(fn, 'wb') as fd:
            fd.write(content)

# vim: set ts=4 sw=4 sts=4 tw=0:

