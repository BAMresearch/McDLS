# -*- coding: utf-8 -*-
# appversion.py

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

# vim: set ts=4 sw=4 sts=4 tw=0:

