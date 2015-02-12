# -*- coding: utf-8 -*-
# gui/liststyle.py

from __future__ import absolute_import # PEP328
import logging

from gui.qt import QtCore, QtGui
from QtCore import QFileInfo
from QtGui import (QWidget, QPalette)

# required for svg graphics support
from gui.qt import QtSvg, QtXml, pluginDirs
from main import makeAbsolutePath

def makeAlternatingRowColorsTransparent(widget):
    palette = widget.palette()
    color = palette.color(QPalette.AlternateBase)
    color.setAlphaF(0.4)
    palette.setColor(QPalette.AlternateBase, color)
    widget.setPalette(palette)

def setBackgroundStyleSheet(widget, imgpath):
    assert isinstance(widget, QWidget)
    makeAlternatingRowColorsTransparent(widget.listWidget)
    stylesheet = """
        #listWidget {{
            background-image:       url({path});
            background-repeat:      no-repeat;
            background-position:    center center;
            background-attachment:  fixed;
            background-color:       white;
        }}
    """
    # convert path to qt style formatting (separators, ...)
    imgpath = QFileInfo(makeAbsolutePath(imgpath)).absoluteFilePath()
    widget.setStyleSheet(stylesheet.format(path = imgpath))

# vim: set ts=4 sts=4 sw=4 tw=0:
