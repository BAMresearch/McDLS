# -*- coding: utf-8 -*-
# gui/version.py

from __future__ import absolute_import # PEP328
from gui.utils.appversion import QAppVersion

version = QAppVersion(
    programName = "McSAS",
    versionNumber = "DLS Spring Edition Mk.VI",
    organizationName = "BAM",
    defaultSettings = dict(
        geometry = "AdnQywABAAAAAAIyAAAA+gAABG8AAANMAAACMwAAARIAAARuAAADSwAAAAAAAA==")
)

