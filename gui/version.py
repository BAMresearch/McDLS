# -*- coding: utf-8 -*-
# gui/version.py

from __future__ import absolute_import # PEP328
from gui.utils.appversion import QAppVersion

version = QAppVersion(
    programName = "McDLS",
    versionNumber = "1.0",
    organizationName = "BAM",
    defaultSettings = dict(
        geometry = "AdnQywABAAAAAAIyAAAA+gAABG8AAANMAAACMwAAARIAAARuAAADSwAAAAAAAA==")
)

