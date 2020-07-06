# -*- coding: utf-8 -*-
# gui/bases/mainwindow/mainwindow.py

from builtins import str
import sys
import logging
from QtWidgets import QMainWindow
from QtCore import QSettings, QByteArray, QTimer
from gui.utils.appversion import QAppVersion
from gui.utils.translate import tr
from gui.utils.signal import Signal
from gui.bases.mainwindow.ui_mainwindow import Ui_MainWindow
from gui.bases.mixins import AppSettings

class MainWindow(QMainWindow, Ui_MainWindow, AppSettings):
    """
    Main window base class.

    Provides functionality for storing and loading application settings,
    managing widgets.
    """
    onStartupSignal = Signal()
    _appversion = None # QAppVersion

    def __init__(self, appversion, parent = None):
        QMainWindow.__init__(self, parent)
        assert QAppVersion.isValid(appversion), "Please provide a valid QAppVersion."
        self.setWindowTitle("{name} {number}"
                .format(name = appversion.name(),
                        number = appversion.number()))
        self.appSettings = QSettings(appversion.organizationName(),
                                     appversion.settingsKey())
        self._appversion = appversion
        self.setupUi(self)
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

    def storeSettings(self):
        if self.appSettings is None:
            return
        geometry = self.saveGeometry()
        windowState = self.saveState()
        # print >>sys.__stderr__, geometry.toBase64()
        # print >>sys.__stderr__, windowState.toBase64()
        self.setRootGroup()
        self.appSettings.setValue("geometry", geometry)
        self.appSettings.setValue("windowState", windowState)

    def restoreSettings(self):
        """
        Load defaults for settings if missing and available.
        """
        if self.appSettings is None:
            return
        defaultSettings = self._appversion.defaultSettings()
        logmsg = "Loaded settings."
        if defaultSettings is not None:
            # loading default settings if not set previously
            custom, default = [], []
            for key, value in defaultSettings.items():
                if self.appSettings.contains(key):
                    custom.append(key)
                else: # qsettings doesn't contain the key, add it
                    self.appSettings.setValue(key,
                            QByteArray.fromBase64(str.encode(value)))
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
                func(self.appSettings.value(key).toByteArray())
            except AttributeError:
                func(self.appSettings.value(key))

# TODO: tests?

# vim: set ts=4 sw=4 sts=4 tw=0:
