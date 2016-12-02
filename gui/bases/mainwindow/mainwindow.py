# -*- coding: utf-8 -*-
# gui/bases/mainwindow/mainwindow.py

from __future__ import absolute_import # PEP328
from builtins import str
import sys
import logging
from QtGui import QMainWindow
from QtCore import QSettings, QByteArray, QTimer
from gui.utils.appversion import QAppVersion
from gui.utils.translate import tr
from gui.utils.signal import Signal
from gui.bases.mainwindow.ui_mainwindow import Ui_MainWindow

class MainWindow(QMainWindow, Ui_MainWindow):
    """
    Main window base class.

    Provides functionality for storing and loading application settings,
    managing widgets.
    """
    onStartupSignal = Signal()
    _appversion = None # QAppVersion
    _appsettings = None # AppSettings

    def __init__(self, appversion, parent = None):
        QMainWindow.__init__(self, parent)
        self.setupUi(self)
        assert QAppVersion.isValid(appversion), "Please provide a valid QAppVersion."
        self._appversion = appversion
        self.setWindowTitle("{name} {number}"
                .format(name = appversion.name(),
                        number = appversion.number()))
        self.restoreSettings()

    def getCommandlineArguments(self):
        """Get command line arguments, if any."""
        arguments = sys.argv[1:]
        programName = self._appversion.name().lower()
        if (len(arguments) > 0 and
            str(arguments[0]).lower().find(programName + ".") >= 0):
            arguments = arguments[1:]
        return arguments

    def show(self):
        QMainWindow.show(self)
        QTimer.singleShot(300, self.onStartupSignal.emit)

    def closeEvent(self, event):
        self.storeSettings()
        QMainWindow.closeEvent(self, event)

    def appSettings(self):
        return self._appsettings

    def storeSettings(self):
        if self._appsettings is None:
            return
        geometry = self.saveGeometry()
        windowState = self.saveState()
        # print >>sys.__stderr__, geometry.toBase64()
        # print >>sys.__stderr__, windowState.toBase64()
        self._appsettings.setValue("geometry", geometry)
        self._appsettings.setValue("windowState", windowState)

    def restoreSettings(self):
        """
        Load defaults for settings if missing and available.
        """
        if self._appsettings is None:
            self._appsettings = QSettings(
                    self._appversion.organizationName(),
                    self._appversion.settingsKey())
        defaultSettings = self._appversion.defaultSettings()
        logmsg = "Loaded settings."
        if defaultSettings is not None:
            # loading default settings if not set previously
            custom, default = [], []
            for key, value in defaultSettings.items():
                if self._appsettings.contains(key):
                    custom.append(key)
                else: # qsettings doesn't contain the key, add it
                    self._appsettings.setValue(key, 
                            QByteArray.fromBase64(value))
                    default.append(key)
            if len(default) > 0:
                logmsg += " Defaults for {0}.".format(", ".join(default))
            if len(custom) > 0:
                logmsg += " Previous ones for {0}.".format(", ".join(custom))
        logging.info(logmsg)
        # restore actual settings
        for key, func in (('geometry', self.restoreGeometry),
                          ('windowState', self.restoreState)):
            try: #QVariant
                func(self._appsettings.value(key).toByteArray())
            except AttributeError:
                func(self._appsettings.value(key))

# TODO: tests?

# vim: set ts=4 sw=4 sts=4 tw=0:
