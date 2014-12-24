# -*- coding: utf-8 -*-
# gui/bases/logwidget.py

from __future__ import absolute_import # PEP328
import os.path
import sys
import re
from gui.qt import QtCore, QtGui
from QtCore import Qt, QRegExp, QTimer
from QtGui import (QApplication, QKeySequence, QTextBrowser, QDesktopServices)
from gui.bases.mixins.titlehandler import TitleHandler
from gui.bases.mixins.contextmenuwidget import ContextMenuWidget
import log
from utils import isString
from utils.lastpath import LastPath
from gui.utils.translate import tr
from gui.utils.filedialog import getSaveFile

# precompile regular expression for urls
URLREGEXP = re.compile(r"((?:[a-z][\w-]+:(?:/{1,3}|[a-z0-9%])|www\d{0,3}[.]|[a-z0-9.\-]+[.][a-z]{2,4}/)(?:[^\s()<>]+|\(([^\s()<>]+|(\([^\s()<>]+\)))*\))+(?:\(([^\s()<>]+|(\([^\s()<>]+\)))*\)|[^\s`!()\[\]{};:'\".,<>?«»“”‘’]))")

def url2href(text):
    """http://daringfireball.net/2010/07/improved_regex_for_matching_urls"""
    return URLREGEXP.sub(r'<a href="\1">\1</a>', text)

class LogWidget(QTextBrowser, ContextMenuWidget):
    """
    Simple TextEdit which can save its contents to file.

    Fill it with text.
    >>> from gui.utils.dialoginteraction import DialogInteraction
    >>> from gui.widgets.logwidget import LogWidget
    >>> testdata = 'blablubba\\nblubb'
    >>> te = DialogInteraction.instance(LogWidget)
    >>> te.append(testdata)
    >>> te.contents() == testdata
    True

    Call save to file dialog and save.
    >>> fn = str(DialogInteraction.query(fileDialogType(), te.saveToFile,
    ...                                  slot = 'accept'))

    Verify written data and remove test file.
    >>> import os
    >>> data = None
    >>> if fn is not None:
    ...     data = open(fn).read()
    ...     os.remove(fn)
    >>> data == testdata
    True

    PlainTextEdit Tests:
    Simple TextEdit which can save its contents to file.

    Fill it with text.
    >>> from utilsgui import DialogInteraction, fileDialogType
    >>> from plaintextedit import PlainTextEdit
    >>> testdata = 'blablubba\\nblubb'
    >>> te = DialogInteraction.instance(PlainTextEdit)
    >>> te.appendPlainText(testdata)
    >>> te.toPlainText() == testdata
    True

    Call save to file dialog and save.
    >>> fn = str(DialogInteraction.query(fileDialogType(), te.saveToFile,
    ...                                  slot = 'accept'))

    Verify written data and remove test file.
    >>> import os
    >>> data = None
    >>> if fn is not None:
    ...     data = open(fn).read()
    ...     os.remove(fn)
    >>> data == testdata
    True
    """
    title = None
    _copyAvailable = None # used for updating context menu accordingly
    _appversion = None
    _buffer = None
    _timer = None
    _updateInterval = 300 # ms
#    _lineBreak = "<br />\n"
    _lineBreak = "\n"

    def __init__(self, parent = None, appversion = None):
        super(LogWidget, self).__init__(parent = parent)
        ContextMenuWidget.__init__(self)
        self.title = TitleHandler.setup(self, "Log")
        self.appversion = appversion
        self.setUndoRedoEnabled(False)
        self.setReadOnly(True)
        self.setTextInteractionFlags(
            Qt.LinksAccessibleByMouse|
            Qt.TextSelectableByMouse)
        try:
            self.setOpenExternalLinks(True)
            self.anchorClicked.connect(QDesktopServices.openUrl)
            self.setOpenLinks(False)
        except: pass
        self._buffer = []
        self._bufferDirty = False
        self._timer = QTimer(self)
        self._timer.setInterval(self._updateInterval)
        self._timer.timeout.connect(self._updateContents)
        log.replaceHandler(log.WidgetHandler(self))
        self._setupActions()
        self._copyAvailable = False
        self.copyAvailable.connect(self.setCopyAvailable)
        self._timer.start()
        QTimer.singleShot(50, self._updateContents)

    @property
    def appversion(self):
        return self._appversion

    @appversion.setter
    def appversion(self, appversion):
        self._appversion = appversion

    def setCopyAvailable(self, yes):
        self._copyAvailable = yes
        self.updateMenu()

    def isCopyAvailable(self):
        return self._copyAvailable

    def _setupActions(self):
        # copy relevant actions from standard context menu
        menu = self.createStandardContextMenu()
        for a in menu.actions():
            title = unicode(a.text()).replace('&', '')
            startswith = title.startswith
            if title.startswith("Copy\t"):
                self.addMenuEntryAction("copy", a, "isCopyAvailable")
            elif title.startswith("Select All\t"):
                self.addMenuEntryAction("selectall", a, "*")
        self.addMenuSeparator()
        self.addMenuEntry("savetofile", tr("Save Log to File"),
                          menuStates = "*", callbacks = self.saveToFile)
        self.updateMenu()

    def isEmpty(self):
        return (len(self.contents()) <= 0)

    def clear(self):
        super(LogWidget, self).clear()
        self._buffer = []


    def _updateContents(self):
        if not self._bufferDirty:
            return
        # remember slider position if not at the bottom
        scroll = self.verticalScrollBar()
        sliderPosition = None
        if scroll.sliderPosition() != scroll.maximum():
            sliderPosition = scroll.sliderPosition()
        # update text
        super(LogWidget, self).append(
                "<pre style=\"white-space: pre-wrap; margin: 0;\">" +
                    self._lineBreak.join(filter(len, self._buffer))
                + "</pre>"
        )
        self._buffer = []
        self._bufferDirty = False
        # restore slider position eventually
        if sliderPosition is not None:
            scroll.setValue(sliderPosition)
        else:
            scroll.setValue(scroll.maximum())

    def scrollToTop(self):
        scroll = self.verticalScrollBar()
        scroll.setValue(scroll.minimum())

    def scrollToBottom(self):
        scroll = self.verticalScrollBar()
        scroll.setValue(scroll.maximum())

    def append(self, text):
        """Appends a new line of text."""
        wasEmpty = self.document().isEmpty()
        self._buffer.extend(url2href(unicode(text.replace('\r', ''))).splitlines())
        self._bufferDirty = True
        if wasEmpty:
            self._updateContents()

    def contents(self):
        return unicode(self.toPlainText())

    def saveToFile(self, filename = None):
        fn = filename
        if not isString(fn):
            filefilter = ["Text files (*.txt)",]
            name, number = "", ""
            if self.appversion:
                name, number = self.appversion.name(), self.appversion.number()
            fnFormat = log.timestamp() + "_{0}-{1}_log.txt"
            fn = fnFormat.format(name, number)
            fn = os.path.join(LastPath.get(), fn)
            fn = getSaveFile(self, "Select a file to save the log to.",
                             fn, filefilter)
        if fn is None or len(fn) < 1:
            return
        with open(fn, 'w') as fd:
            fd.write(bytearray(self.contents(), "utf8"))
        return fn

    def onCloseSlot(self):
        """Workaround to exit a program during calculation.
        Logging a message processes Qt events (see append()).
        It makes sure this slot is called if connected to a application close
        signal. For example, emitted by the main window on closeEvent()."""
        QApplication.instance().exit()

# on single file execution, run doctests
if __name__ == "__main__":
    import doctest
    doctest.testmod()

# vim: set ts=4 sts=4 sw=4 tw=0:
