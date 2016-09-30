# -*- coding: utf-8 -*-
# gui/version.py

from __future__ import absolute_import # PEP328
from gui.utils.appversion import QAppVersion

version = QAppVersion(
    programName = "McSAS",
    versionNumber = "DLS Edition Mk.VIII",
    organizationName = "BAM",
    defaultSettings = dict(
        geometry = "AdnQywABAAAAAAIyAAAA+gAABG8AAANMAAACMwAAARIAAARuAAADSwAAAAAAAA==")
)

