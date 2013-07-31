# -*- coding: utf-8 -*-
# logwidget.py

import os.path
from cutesnake.qt import QtCore, QtGui
from QtCore import Qt
from QtGui import QKeySequence, QTextEdit
from cutesnake.widgets.mixins.titlehandler import TitleHandler
from cutesnake.widgets.mixins.contextmenuwidget import ContextMenuWidget
from cutesnake.log import Log, WidgetHandler
from cutesnake.utils import isString
from cutesnake.utils.lastpath import LastPath
from cutesnake.utils.translate import tr
from cutesnake.utilsgui.filedialog import getSaveFile
from cutesnake.utilsgui import processEventLoop

class LogWidget(QTextEdit, ContextMenuWidget):
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
    _copyAvailable = None # used for updating context menu accordingly
    _appversion = None
    title = None

    def __init__(self, parent = None, appversion = None):
        QTextEdit.__init__(self, parent)
        ContextMenuWidget.__init__(self)
        self.title = TitleHandler.setup(self, "Log")
        self.appversion = appversion
        self.setUndoRedoEnabled(False)
        self.setReadOnly(True)
        self.setTextInteractionFlags(
            Qt.LinksAccessibleByMouse|
            Qt.TextSelectableByMouse)
        Log.setHandler(WidgetHandler(self))
        self._setupActions()
        self._copyAvailable = False
        self.copyAvailable.connect(self.setCopyAvailable)

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
            title = a.text().remove('&')
            if title.startsWith("Copy\t"):
                self.addMenuEntryAction("copy", a, "isCopyAvailable")
            elif title.startsWith("Select All\t"):
                self.addMenuEntryAction("selectall", a, "*")
        self.addMenuSeparator()
        self.addMenuEntry("savetofile", tr("Save Log to File"),
                          menuStates = "*", callbacks = self.saveToFile)
        self.updateMenu()

    def isEmpty(self):
        return (len(self.contents()) <= 0)

    def append(self, text):
        """Interface method for WidgetHandler."""

        if text[-1] == '\n':
            text = text[:-1]

        html = self.toHtml()
        if self.contents().isEmpty():
            html = ""
        if getattr(self, "hadCR", False):
            # handle carriage return
            lastbreak = html.lastIndexOf('\n')
            if lastbreak >= 0:
                html.truncate(lastbreak)
            scroll = self.verticalScrollBar()
            scroll.setValue(scroll.maximum())
            self.hadCR = False

        if text[-1] == '\r':
            self.hadCR = True

        super(LogWidget, self).setHtml(
                html + text.replace('\r', '').replace('\n','<br />\n')
        )
        processEventLoop()

    def contents(self):
        return self.toPlainText()

    def saveToFile(self, filename = None):
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
            fd.write(self.contents())
        return fn

# vim: set ts=4 sts=4 sw=4 tw=0:
