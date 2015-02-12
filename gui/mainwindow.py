# -*- coding: utf-8 -*-
# gui/mainwindow.py

from __future__ import absolute_import # PEP328
import os.path
import re
import sys
import logging

from gui.qt import QtCore, QtGui
from gui.utils.signal import Signal
from QtCore import Qt, QSettings, QFileInfo
from QtGui import (QWidget, QHBoxLayout, QVBoxLayout, QPushButton,
                   QSizePolicy, QApplication, QToolBox, QIcon)
from gui.bases.mainwindow import MainWindow as MainWindowBase
from gui.bases.logwidget import LogWidget
from gui.bases.dockwidget import DockWidget
from gui.utils.filedialog import getOpenFiles
from utils.lastpath import LastPath
from utils import isMac
from gui.utils import processEventLoop
from gui.utils.displayexception import DisplayException
from gui.version import version
from gui.calc import Calculator
from sasdata import SASData

from gui.scientrybox import SciEntryBox

INFOTEXT = """One or more selected files are read in and passed to Brian Pauws Monte-Carlo size distribution analysis program for 1D SAXS data.

The convergence criterion can be set by the user. If it is not reached no output is generated, just this log is saved to file. On success, the resulting size distribution and the data fit are stored to files with uncertainties.

Output files start with the base name of the input file. They have the current date+time appended to avoid overwriting existing results.

[ For more information, please see http://www.mcsas.net ]"""

CHANGESTEXT = (u"""
Latest changes:
- Number input boxes now allow scientific notation input
- Stability improvements and code cleanup
- Improved plotting routine
- Range statistics list functionality added
- Improved division of GUI items into vertical tabs
- More information shown in data list
- Correct handling when fitting multiple files
- Range estimate now uses minimum Q by itself, no longer considers Q-spacing
- LMA Dense Spheres, spherical and ellipsoidal core shell models work again
- Extended internal parameter functionality, using JSON defaults file
- Improvements towards implementing RangeInfo in the GUI

Changes in 0.0.11:
- distribution statistics log output and writing to stats file
- plain SAS evaluation (no fit) if no param is active
- all output files with extensions and storing in settings
- configuration file now *.cfg
- single start/stop button

Changes in 0.0.10:
- data file listing widget on top with short info
  - double-click uses sphere size est. for sphere model radius range
- start/stop buttons
  - allows to abort current calculation and restart
- column names in output files
- data file names stored with settings
- switch to enable/disable background level fit
- multiple plot figures on Windows supported

Changes in 0.0.9:
- added GUI to public McSAS repository

Changes in 0.0.8:
- new model: GaussianChain, verified against SASfit:
  http://sasfit.sf.net/manual/Gaussian_Chain#Gauss_2
- fixed volume function exponent in Kholodenko
  was v² instead of just v

Changes in 0.0.7:
- Using restructured McSAS
- building GUI for models and global settings dynamically
- new model: Kholodenkos worm-like structure, verified against SASfit

Changes in 0.0.6:
Fixed handling of negative values in PDH data files.

Changes in 0.0.5:
'Number-weighted distributions now come with correct-looking observability limits.'
 https://bitbucket.org/pkwasniew/mcsas/commits/81bbf84

""".replace('\n\n', '<hr />'))
CHANGESTEXT = re.sub(r"([\s\w]*[cC]hanges.*\:)",
                     r"<strong>\1</strong>",
                     CHANGESTEXT)

# required for svg graphics support
from gui.qt import QtSvg, QtXml, pluginDirs
from gui.rangelist import RangeList
from gui.settingswidget import SettingsWidget
from gui.algorithmwidget import AlgorithmWidget
from gui.modelwidget import ModelWidget
from gui.filelist import FileList
from main import makeAbsolutePath

def eventLoop(args):
    """Starts the UI event loop and get command line parser arguments."""
    app = QApplication(sys.argv)
    for pluginDir in pluginDirs(): # required for svg graphics support
        app.addLibraryPath(pluginDir)
    mw = MainWindow(args = args)
    mw.show()
    return app.exec_()


class ToolBox(QToolBox):
    """QToolBox containing the widgets for user settings.
    Used to propagate resize events to child widgets to enable responsive behaviour.
    On MacOS, fixes failed detection of size changes in child widget due to scroll area.
    """
    def resizeEvent(self, event):
        child = self.currentWidget()
        if isinstance(child, AlgorithmWidget):
            child.resizeEvent(event)

class MainWindow(MainWindowBase):
    onCloseSignal = Signal()
    _args = None # python command line arguments parser
    _calculator = None # calculator calling algorithm on all data

    def __init__(self, parent = None, args = None):
        # calls setupUi() and restoreSettings()
        MainWindowBase.__init__(self, version, parent)
        self._args = args

    @property
    def calculator(self):
        """Returns a calculator object."""
        if self._calculator is None:
            self._calculator = Calculator()
        return self._calculator

    def setupUi(self, *args):
        # called in MainWindowBase.__init__()
        # put the log widget at the bottom
        self.addDockWidget(Qt.BottomDockWidgetArea, self._setupLogWidget())
        # file widget at the top
        self.toolbox = ToolBox(self)
        self._addToolboxItem(self._setupFileWidget())
        self._addToolboxItem(self._setupAlgoWidget())
        self._addToolboxItem(self._setupModelWidget())
        self._addToolboxItem(self._setupStatsWidget())

        # set up central widget of the main window
        self.centralLayout = QVBoxLayout()
        # put buttons in central widget
        self.centralLayout.addWidget(self.toolbox)
        self.centralLayout.addWidget(self._setupStartButton())
        centralWidget = QWidget(self)
        centralWidget.setLayout(self.centralLayout)
        centralWidget.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        self.setCentralWidget(centralWidget)
        self.onStartupSignal.connect(self.initUI)
        # set program icon, same for Win+Lin
        icopath = "resources/icon/mcsas.ico"
        if isMac():
            icopath = "resources/icon/mcsas.icns"
        icopath = QFileInfo(makeAbsolutePath(icopath)).absoluteFilePath()
        self.setWindowIcon(QIcon(icopath))

    def _addToolboxItem(self, widget):
        self.toolbox.addItem(widget, "{n}. {t}"
                                     .format(n = self.toolbox.count()+1,
                                             t = widget.title()))
        try:
            widget.layout().setContentsMargins(0, 0, 0, 0)
        except: pass

    def _setupFileWidget(self):
        # set up file widget
        fileWidget = FileList(self, title = "Data Files",
                              withBtn = False, nestedItems = False)
        fileWidget.setHeader(SASData.displayDataDescr)
        fileWidget.setToolTip(
                "Right-click to add datafiles.\n" +
                "Double click to use the estimated size for the model.")
        self.fileWidget = fileWidget
        return fileWidget

    def _setupAlgoWidget(self):
        """Set up property widget with settings."""
        self.algoWidget = AlgorithmWidget(self, self.calculator)
        return self.algoWidget

    def _setupModelWidget(self):
        """Set up property widget with settings."""
        self.modelWidget = ModelWidget(self, self.calculator)
        self.fileWidget.sigSphericalSizeRange.connect(
                self._onSphericalSizeRange)
        return self.modelWidget

    def _onSphericalSizeRange(self, *args):
        self.toolbox.setCurrentWidget(self.modelWidget)
        self.modelWidget.setSphericalSizeRange(*args)

    def _setupStatsWidget(self):
        """Set up property widget with settings."""
        # setup similar to the file widget
        self.statsWidget = RangeList(parent = self,
                                     calculator = self.calculator,
                                     title = "Post-fit Analysis",
                                     withBtn = False, nestedItems = False)
        self.modelWidget.setStatsWidget(self.statsWidget)
        self.modelWidget.sigModelChanged.connect(self.statsWidget.updateHistograms)
        return self.statsWidget

    def _setupLogWidget(self):
        """Set up widget for logging output."""
        logDock = DockWidget(self, LogWidget, appversion = version)
        logWidget = logDock.child
        self.onCloseSignal.connect(logWidget.onCloseSlot)
        logWidget.setSizePolicy(QSizePolicy.Preferred,
                                QSizePolicy.Expanding)
        logWidget.append(INFOTEXT)
        if len(CHANGESTEXT):
            logWidget.append(CHANGESTEXT)
        logWidget.append("\n\r")
        self.logWidget = logWidget
        return logDock

    def _setupStartButton(self):
        """Set up "Start/Stop" - button."""
        self.startStopBtn = QPushButton()
        self.startStopBtn.setCheckable(True)
        self.startStopBtn.clicked[bool].connect(self.onStartStopClick)
        btnLayout = QHBoxLayout()
        btnLayout.setContentsMargins(0, 0, 0, 0)
        self.startStopBtn.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Maximum)
        btnLayout.addWidget(self.startStopBtn)
        btnWidget = QWidget(self)
        btnWidget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Maximum)
        btnWidget.setLayout(btnLayout)
        return btnWidget

    def restoreSettings(self):
        MainWindowBase.restoreSettings(self)
        settings = self.appSettings()
        for settingsWidget in self.algoWidget, self.modelWidget:
            settingsWidget.appSettings = self.appSettings()
            settingsWidget.restoreSession()
        try:
            value = unicode(settings.value("lastpath").toString())
        except AttributeError: # QVariant
            value = unicode(settings.value("lastpath"))
        if os.path.isdir(value):
            LastPath.set(value)

    def storeSettings(self):
        MainWindowBase.storeSettings(self)
        settings = self.appSettings()
        for settingsWidget in self.algoWidget, self.modelWidget:
            settingsWidget.storeSession()
        settings.setValue("lastpath", LastPath.get())
        settings.sync()
        return
        # test for additionally storing settings to file
        tempSettings = QSettings("/tmp/qsettings.test", QSettings.IniFormat)
        for key in settings.allKeys():
            if key in ('geometry', 'windowState', 'lastpath'):
                continue
            tempSettings.setValue(key, settings.value(key))
        tempSettings.sync()

    def fileDialog(self):
        filenames = getOpenFiles(self, "Load one or more data files",
                                 LastPath.get(), multiple = True)
        self.loadFiles(filenames)

    def initUI(self):
        self.logWidget.scrollToTop()
        self.fileWidget.loadData(getattr(self._args, "fnames", []))
        self.onStartStopClick(getattr(self._args, "start", False))

    def _updateWidgets(self):
        for w in self.findChildren(SettingsWidget):
            w.updateAll()

    def onStartStopClick(self, checked):
        processEventLoop()
        if checked:
            self.startStopBtn.setText("stop")
            self.startStopBtn.setChecked(True)
            self._updateWidgets() # get latest input in case sth didn't update
            self.calc()
        # run this also for 'start' after calculation
        self.calculator.stop()
        self.startStopBtn.setText("start")
        self.startStopBtn.setChecked(False)

    def calc(self):
        if len(self.fileWidget) <= 0:
            return
        self.logWidget.clear()
        self.logWidget.scrollToBottom()
        self.fileWidget.updateData(updateFunc = self.calculator,
                                   stopFunc = self.calculator.isStopped,
                                   showProgress = False)

    def closeEvent(self, closeEvent):
        super(MainWindow, self).closeEvent(closeEvent)
        self.onStartStopClick(False)
        self.onCloseSignal.emit()

    def keyPressEvent(self, keyEvent):
        if keyEvent.key() == Qt.Key_Escape and self.startStopBtn.isChecked():
            self.onStartStopClick(False) # hit 'stop'
            logging.info("Calculation aborted by user interupt!")

# vim: set ts=4 sts=4 sw=4 tw=0:
