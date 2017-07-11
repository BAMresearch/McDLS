# -*- coding: utf-8 -*-
# gui/mainwindow.py

from __future__ import absolute_import # PEP328
from builtins import str
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
from dataobj import SASData

from gui.scientrybox import SciEntryBox

# required for svg graphics support
from gui.qt import QtSvg, QtXml, pluginDirs
from gui.rangelist import RangeList
from gui.algorithmwidget import AlgorithmWidget
from gui.datawidget import DataWidget
from gui.optimizationwidget import OptimizationWidget
from gui.modelwidget import ModelWidget
from gui.filelist import FileList
from main import makeAbsolutePath

INFOTEXT = (u"""
{programName}:\
 Software for the retrieval of model parameter distributions from\
 dynamic light scattering data.""".format(programName = version.name())
+ u"""
Output files of a Monte Carlo run are stored in a directory named after\
 the input file followed with a timestamp to avoid overwriting existing\
 results.

<br /><strong>Literature:</strong>
This suite and method are detailed in:
- Bressler, I. et al, J. Appl. Cryst. 00: 000-111.
  http://dx.doi.org/
- Bressler, I. et al, J. Appl. Cryst. 48: 962-969.
  http://dx.doi.org/10.1107/S1600576715007347
- Pauw, B. R. et al., J. Appl. Cryst. 46: 365-371.
  http://dx.doi.org/10.1107/S0021889813001295

<br /><strong>What to do in case of unsuccessful fits:</strong>
If convergence is not reached no output is generated, and only this
log is stored in a file. On success, the resulting size distribution
and data fit are stored to files with uncertainties.

If convergence is not reached, the following can be attempted:

1) Adjust the convergence criterion to a larger value. 
   The convergence criterion can be adjusted by the user, to support
   data whose uncertainties are too large or too small.
2) Verify that the parameter range of the model is appropriate. Very
   wide ranges may prevent the correct solution from appearing within
   the limited number of iterations attempted by the program.

<style>li { margin: .5em; }</style>
<br /><strong>Notes on DLS data handling:</strong><ul>
<li>The program expects *.ASC data files, containing "ALV-7004 CGS-8F Data"\
 in the first line of the file. It shows a warning otherwise.</li>

<li>It is recommended to load multiple files from measurements of the same\
 sample at once. They are averaged automatically in order to infer\
 uncertainties for all values. In a second step, the averaged data set is\
 split into one data set for every angle. For example, loading six files of\
 measurements at eight angles will result in eight averaged data sets being\
 listed in the program.</li>

<li>Files from more than one sample can be loaded, all files sharing\
 the same sample name are averaged along matching angles.</li>

<li>The measurement numbers ("indices") are taken from the file name of\
 each data file. They consist of a group index and a measurement index which\
 are shown in the column "Measurements".</li>

<li>The program processes and fits only those data sets which are selected.\
 If none are selected, it processes all data sets. The selection of all data\
 sets can be toggled by <Ctrl-A> key presses or via context menu at the \
 data set list (right mouse button).</li>

<li>The menu "Data Settings" shows configuration options for the first of\
 the selected data sets. Changes are transferred to all data sets sharing its\
 sample name.</li>

<li>The option "Calc. series statistics" in the optimization menu combines the\
 distribution moments of all data sets in one output file\
 ("[...]_seriesStats.dat") for each active parameter. Finally, it shows a\
 simple plot of the mean and its standard deviation across\
 the scattering angles.</li>
</ul>
[ For more information, please see http://www.mcdls.net ]
""")

CHANGESTEXT = (u"""

""".replace('\n\n', '<hr />')) # html horz. line instead of 2x newline
# make 'changes xyz' bold; wondering how much weight markdown might add (?)
CHANGESTEXT = re.sub(r"([\s\w]*[cC]hanges.*\:)",
                     r"<strong>\1</strong>",
                     CHANGESTEXT)

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
        if isinstance(child, OptimizationWidget):
            child.resizeEvent(event)
        if isinstance(child, DataWidget):
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
        self._addToolboxItem(self._setupDataWidget())
        self._addToolboxItem(self._setupOptimWidget())
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
        self.onStartupSignal.connect(self.initUi)
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
        fileWidget.setToolTip(
                "Right-click to add datafiles.\n" +
                "Double click to use the estimated size for the model.")
        self.fileWidget = fileWidget
        return fileWidget

    def _setupDataWidget(self):
        """Set up property widget with settings."""
        self.dataWidget = DataWidget(self, self.appSettings)
        self.dataWidget.sigConfig.connect(self.fileWidget.setDataConfig)
        self.fileWidget.sigSelectedData.connect(self.dataWidget.onDataSelected)
        self.fileWidget.sigEmpty.connect(self.dataWidget.onEmptyDataList)
        return self.dataWidget

    def _setupOptimWidget(self):
        """Set up property widget with settings."""
        self.optimWidget = OptimizationWidget(self, self.calculator.algo, self.appSettings)
        self.fileWidget.sigSelectedData.connect(self.optimWidget.onDataSelected)
        return self.optimWidget

    def _setupModelWidget(self):
        """Set up property widget with settings."""
        self.modelWidget = ModelWidget(self, self.calculator, self.appSettings)
        self.fileWidget.sigSphericalSizeRange.connect(
                self._onSphericalSizeRange)
        self.fileWidget.sigSelectedData.connect(self.modelWidget.onDataSelected)
        return self.modelWidget

    def _onSphericalSizeRange(self, *args):
        if self.modelWidget.setSphericalSizeRange(*args):
            self.toolbox.setCurrentWidget(self.modelWidget)

    def _setupStatsWidget(self):
        """Set up property widget with settings."""
        # setup similar to the file widget
        self.statsWidget = RangeList(parent = self,
                                     calculator = self.calculator,
                                     appSettings = self.appSettings,
                                     title = "Post-fit Analysis",
                                     withBtn = False, nestedItems = False)
        self.modelWidget.setStatsWidget(self.statsWidget)
        self.modelWidget.sigBackendUpdated.connect(self.statsWidget.updateHistograms)
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
        for settingsWidget in (self.optimWidget, self.modelWidget):
            # no self.dataWidget, it is restored on demand internally
            settingsWidget.restoreSession()
        if self.appSettings is None:
            return
        try:
            value = str(self.appSettings.value("lastpath").toString())
        except AttributeError: # QVariant
            value = str(self.appSettings.value("lastpath"))
        if os.path.isdir(value):
            LastPath.set(value)

    def storeSettings(self):
        MainWindowBase.storeSettings(self)
        for settingsWidget in (self.optimWidget, self.modelWidget,
                               self.dataWidget):
            settingsWidget.storeSession()
        if self.appSettings is not None:
            self.appSettings.setValue("lastpath", LastPath.get())
            self.appSettings.sync()
        return
        # test for additionally storing settings to file
        tempSettings = QSettings("/tmp/qsettings.test", QSettings.IniFormat)
        for key in self.appSettings.allKeys():
            if key in ('geometry', 'windowState', 'lastpath'):
                continue
            tempSettings.setValue(key, self.appSettings.value(key))
        tempSettings.sync()

    def fileDialog(self):
        filenames = getOpenFiles(self, "Load one or more data files",
                                 LastPath.get(), multiple = True)
        self.loadFiles(filenames)

    def initUi(self):
        self.logWidget.scrollToTop()
        self.fileWidget.loadData(getattr(self._args, "fnames", []))
        self.onStartStopClick(getattr(self._args, "start", False))

    def _updateWidgetsFinally(self):
        for w in self.findChildren(AlgorithmWidget):
            w.updateAll()
        if (not len(self.statsWidget.data())
            and len(self.modelWidget.model.activeParams())):
            # make sure there is an histogram range defined,
            # otherwise ask the user for one
            self.toolbox.setCurrentWidget(self.statsWidget)
            self.statsWidget.loadData()

    def onStartStopClick(self, checked):
        processEventLoop()
        if checked:
            # # write HDF datafile
            # self.hdfStore("test3.h5")
            self.startStopBtn.setText("stop")
            self.startStopBtn.setChecked(True)
            self._updateWidgetsFinally() # get latest input in case sth didn't update
            self.calc()
        # run this also for 'start' after calculation
        self.calculator.stop()
        self.startStopBtn.setText("start")
        self.startStopBtn.setChecked(False)

    # def hdfWrite(self, hdf):
    #     hdf.writeMember(self, "calculator")

    def calc(self):
        if len(self.fileWidget) <= 0:
            return
        self.logWidget.clear()
        self.logWidget.scrollToBottom()
        idx, data = self.fileWidget.currentSelection()
        # process the selected data only if there is a selection
        selectedOnly = data is not None
        self.calculator.prepare()
        self.fileWidget.updateData(updateFunc = self.calculator,
                                   stopFunc = self.calculator.isStopped,
                                   selectedOnly = selectedOnly,
                                   showProgress = False)
        self.calculator.postProcess()

    def closeEvent(self, closeEvent):
        super(MainWindow, self).closeEvent(closeEvent)
        self.onStartStopClick(False)
        self.onCloseSignal.emit()

    def keyPressEvent(self, keyEvent):
        if keyEvent.key() == Qt.Key_Escape and self.startStopBtn.isChecked():
            self.onStartStopClick(False) # hit 'stop'
            logging.info("Calculation aborted by user interupt!")

# vim: set ts=4 sts=4 sw=4 tw=0:
