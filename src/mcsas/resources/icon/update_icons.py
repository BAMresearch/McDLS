#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# update_icons.py

"""
Updates all icons based on the given SVG file.
Call it like that from the McSAS source tree:

    ./resources/icon/update_icons.py resources/icon/mcsas.svg

"""

import sys
import os.path
from create_icons import createIcons
from build_iconset import buildIconSet

if __name__ == "__main__":
    assert len(sys.argv) > 1, (
        "Please provide an icon SVG file!")
    fn = sys.argv[1]
    createIcons(fn)
    fn = os.path.splitext(fn)[0]
    buildIconSet(fn + ".ico")
    buildIconSet(fn + ".icns")

# vim: set ts=4 sw=4 sts=4 tw=0:
