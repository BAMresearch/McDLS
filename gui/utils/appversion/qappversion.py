# -*- coding: utf-8 -*-
# gui/utils/appversion/qappversion.py

from __future__ import absolute_import # PEP328
from gui.utils.appversion import AppVersion
from utils import isString
from gui.qt import QtCore
from QtCore import QCoreApplication as coreApp

class QAppVersion(AppVersion):
    """
    Set QCoreApplication properties based on version meta data.
    """
    def __init__(self, *args, **kwargs):
        AppVersion.__init__(self, *args, **kwargs)
        self._setApplicationMetaData()

    def settingsKey(self):
        """
        Version dependent settings key.
        """
        majorMinor = ".".join(str(self.number()).split(".")[0:2])
        return "{0}_{1}".format(self.name(), majorMinor)

    def _setApplicationMetaData(self):
        for func, data in (
                (coreApp.setApplicationName, self.name()),
                (coreApp.setApplicationVersion, self.number()),
                (coreApp.setOrganizationName, self.organizationName()),
                (coreApp.setOrganizationDomain, self.organizationDomain())
                ):
            if isString(data):
                func(data)

# vim: set ts=4 sw=4 sts=4 tw=0:

