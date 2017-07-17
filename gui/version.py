# -*- coding: utf-8 -*-
# gui/version.py

from __future__ import absolute_import # PEP328
from gui.utils.appversion import QAppVersion

version = QAppVersion(
    programName = "McDLS",
    versionNumber = "1.0.2",
    organizationName = "BAM",
    colors = "#ff470e;#7c260b;#e55528;#ff6a3d", # inner, outer, circles, shell
    defaultSettings = dict(
        geometry = "AdnQywABAAAAAAIyAAAA+gAABG8AAANMAAACMwAAARIAAARuAAADSwAAAAAAAA==")
)

