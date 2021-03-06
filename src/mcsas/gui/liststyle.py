# -*- coding: utf-8 -*-
# gui/liststyle.py

import logging
from QtCore import QFileInfo
from QtGui import QPalette
from QtWidgets import QWidget

# required for svg graphics support
from .qt import QtSvg, QtXml, pluginDirs
from ..main import makeAbsolutePath

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
            background-color:       {bgcol};
        }}
    """
    # convert path to qt style formatting (separators, ...)
    imgpath = QFileInfo(makeAbsolutePath(imgpath)).absoluteFilePath()
    # also set background color of image to same used in GUI (covers dark mode)
    bgcol = widget.listWidget.palette().color(QPalette.Base)
    widget.setStyleSheet(stylesheet.format(path=imgpath, bgcol=bgcol.name()))

# vim: set ts=4 sts=4 sw=4 tw=0:
