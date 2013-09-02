# -*- coding: utf-8 -*-
# logwidget.py

import os.path
import sys
import re
from cutesnake.qt import QtCore, QtGui
from QtCore import Qt, QRegExp, QTimer
from QtGui import (QApplication, QKeySequence, QTextBrowser, QDesktopServices)
from cutesnake.widgets.mixins.titlehandler import TitleHandler
from cutesnake.widgets.mixins.contextmenuwidget import ContextMenuWidget
from cutesnake.log import Log, WidgetHandler
from cutesnake.utils import isString
from cutesnake.utils.lastpath import LastPath
from cutesnake.utils.translate import tr
from cutesnake.utilsgui.filedialog import getSaveFile
from cutesnake.utilsgui import processEventLoop

import cProfile

# precompile regular expression for urls
URLREGEXP = re.compile(r"((?:[a-z][\w-]+:(?:/{1,3}|[a-z0-9%])|www\d{0,3}[.]|[a-z0-9.\-]+[.][a-z]{2,4}/)(?:[^\s()<>]+|\(([^\s()<>]+|(\([^\s()<>]+\)))*\))+(?:\(([^\s()<>]+|(\([^\s()<>]+\)))*\)|[^\s`!()\[\]{};:'\".,<>?«»“”‘’]))")

def url2href(text):
    """http://daringfireball.net/2010/07/improved_regex_for_matching_urls"""
    return URLREGEXP.sub(r'<a href="\1">\1</a>', text)

class LogWidget(QTextBrowser, ContextMenuWidget):
    """
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
    _updateInterval = 1000 # ms
    _lineBreak = "<br />\n"

    def __init__(self, parent = None, appversion = None):
        QTextBrowser.__init__(self, parent)
        ContextMenuWidget.__init__(self)
        self.title = TitleHandler.setup(self, "Log")
        self.appversion = appversion
        self.setUndoRedoEnabled(False)
        self.setReadOnly(True)
        self.setOpenExternalLinks(True)
        self.setTextInteractionFlags(
            Qt.LinksAccessibleByMouse|
            Qt.TextSelectableByMouse)
        self._buffer = []
        self._bufferDirty = False
        self._timer = QTimer(self)
        self._timer.setInterval(self._updateInterval)
        self._timer.timeout.connect(self._updateContents)
        Log.setHandler(WidgetHandler(self))
        self._setupActions()
        self._copyAvailable = False
        self.copyAvailable.connect(self.setCopyAvailable)
        self.anchorClicked.connect(QDesktopServices.openUrl)
        self.setOpenLinks(False)
        self._timer.start()
        QTimer.singleShot(50, self._updateContents)

        self.profile = cProfile.Profile()

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

    def _appendBuffer(self, text):
        # Handling carriage return \r correctly for overwriting a single line
        # makes this a bit complex ...
        if getattr(self, "lastr", -1) >= 0:
            # delete \r preceding text in the previous line until previous \n
            prev = self._buffer.pop(-1)
            lastn = max(0, prev.rfind('\n', 0, self.lastr))
            prev = prev[lastn-1:self.lastr]
            if len(prev):
                self._buffer.append(prev)
                self._bufferDirty = True
            self.lastr = -1 # reset
        # look for \r in current line
        self.lastr = text.rfind('\r')
        lastn = -1
        if self.lastr >= 0:
            lastn = max(0, text.rfind('\n', 0, self.lastr))
            if len(text[self.lastr+1:]):           # there is text following
                text = text[lastn-1:self.lastr] # delete preceding text
                self.lastr = -1                    # reset
        if text[-1] == '\n':                       # remove trailing newlines
            text = text[:-1]
        self._buffer.append(url2href(text.replace('\n', self._lineBreak)))
        self._bufferDirty = True

    def _updateContents(self):
        logfile("_updateContents")
        if not self._bufferDirty:
            return
        # remember slider position if not at the bottom
        scroll = self.verticalScrollBar()
        sliderPosition = None
        if scroll.sliderPosition() != scroll.maximum():
            sliderPosition = scroll.sliderPosition()
        # update text
        super(LogWidget, self).setHtml(
                self._lineBreak.join(self._buffer)
        )
        self._bufferDirty = False
        # restore slider position eventually
        if sliderPosition is not None:
            scroll.setValue(sliderPosition)
        else:
            scroll.setValue(scroll.maximum())

    def append(self, text):
        """Appends a new line of text."""
        self.profile.enable()
        logfile("append {0}".format(len(self._buffer)))

        self._appendBuffer(unicode(text))
        # process qt events always
        processEventLoop()

        self.profile.disable()

    def contents(self):
        return unicode(self.toPlainText())

    def saveToFile(self, filename = None):
        print "storing profile"
        self.profile.dump_stats("/tmp/profile.log")

        fn = filename
        if not isString(fn):
            filefilter = ["Text files (*.txt)",]
            name, number = "", ""
            if self.appversion:
                name, number = self.appversion.name(), self.appversion.number()
            fnFormat = Log.timestamp() + "_{0}-{1}_log.txt"
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

# vim: set ts=4 sts=4 sw=4 tw=0:
