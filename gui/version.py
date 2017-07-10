# -*- coding: utf-8 -*-
# gui/version.py

from __future__ import absolute_import # PEP328
from gui.utils.appversion import QAppVersion

version = QAppVersion(
    programName = "McDLS",
    versionNumber = "1.0",
    organizationName = "BAM",
    colors = "#29235c;#1d70b7;#1d70b7;#35a8e0", # inner, outer, circles, shell
    defaultSettings = dict(
        geometry = "AdnQywABAAAAAAIyAAAA+gAABG8AAANMAAACMwAAARIAAARuAAADSwAAAAAAAA==")
)

