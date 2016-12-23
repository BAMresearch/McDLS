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
<strong>major changes in McDLS Mk.IX:</strong><ul>
<style>li { margin: .5em; }</style>
<li>alternative Goodness-of-Fit indicator output added</li>
</ul>
<strong>major changes in McDLS Mk.VIII:</strong><ul>
<style>li { margin: .5em; }</style>
<li>fixed histogramming regarding calculated DLS amplitude</li>
<li>option to square the amplitude or not (V²Phi² vs. V·Phi)</li>
<li>option to fix the first point of the model to the measured data</li>
<li>surface weighted histogram</li>
<li>current repetition included in progress log message</li>
<li>Windows: URL links in log window should work again</li>
<li>mostly fixed overlapping histogram plots</li>
<li>fit configuration written to HDF5</li>
<li>some bug fixes</li>
</ul>
<strong>major changes in McDLS Mk.VII:</strong><ul>
<style>li { margin: .5em; }</style>
<li>fix for long path names on Windows, supposed to work with >260 characters now</li>
<li>fixed one possible source for crashes in the UI</li>
<ul><li><span style="text-decoration: underline;">Feel free to report any crashes!</span> At best with a short description of the last action done with the program:</li>
<li>What triggered a crash?</li>
<li>What was clicked last?</li>
<li>Which program version on which Operating System?</li>
</ul>
</ul>
<strong>major changes in McDLS Mk.VI:</strong><ul>
<style>li { margin: .5em; }</style>
<li>volume squared distribution, aka. intensity weighted</li>
<li>fixes for loading one-angled DLS .ASC data, simulated by CNTb/CONTIN</li>
</ul>
<strong>major changes in McDLS Mk.V:</strong><ul>
<style>li { margin: .5em; }</style>
<li>stores also the averaged DLS data along with fit results</li>
<li>scattering angle dependent form factor in DLSSphere model enabled by default</li>
</ul>
<strong>major changes in McDLS Mk.IV:</strong><ul>
<style>li { margin: .5em; }</style>
<li>fixed volume fraction calculus in histogramming for DLS</li>
<li>weighting compensation has to be == 1! (still adjustable for testing)</li>
<li>fixed scattering angle dependent form factor in DLSSphere model</li>
<li>Skips the first 2 data points by default</li>
<li>more informative plot labels in series analysis</li>
</ul>

<strong>Notes on DLS data handling:</strong><ul>
<style>li { margin: .5em; }</style>
<li>The program expects *.ASC data files, containing "ALV-7004 CGS-8F Data" in the first line of the file. It shows a warning otherwise.</li>

<li>It is recommended to load multiple files from measurements of the same sample at once. They are averaged automatically in order to infer uncertainties for all values. In a second step, the averaged data set is separated in one data set for every angle. For example, loading six files of measurements at eight angles will result in eight averaged data sets being listed in the program.</li>

<li>Files from more than one sample can be loaded, all files containing the same sample name and angles are averaged.</li>

<li>The measurement numbers ("indices") are taken from the file name of each data file. They consist of a group index and a measurement index which are shown in the column "Measurements".</li>

<li>The behaviour for selected data sets changed: The program processes and fits only those data sets which are selected. If none are selected it processes all data sets (old behaviour). The selection of all data sets can be toggled by <Ctrl-A> key presses or via context menu at the data set list (right mouse button). In previous versions, the selection of data sets did not change anything.</li>

<li>The new menu "Data Settings" shows configuration options for the first of the selected data sets. Changes are transferred to all selected data sets of the same type (DLS or SAS).</li>

<li>In the optimization menu there is a new option "combine statistics [...]". If enabled, it combines the distribution moments of all data sets in one output file ("[...]_seriesStats.dat") for each active parameter. Finally, it shows a simple plot of the mean and its standard deviation for all angles.</li>

</ul><hr />
McSAS: Software for the retrieval of model parameter distributions from scattering patterns.

Output files of a Monte Carlo run are stored in a directory named after the input file followed with a timestamp to avoid overwriting existing results.

<br /><strong>Literature:</strong>
This suite and method are detailed in:
- Bressler, I. et al, J. Appl. Cryst. 48: 962-969.
  http://dx.doi.org/10.1107/S1600576715007347
- Pauw, B. R. et al., J. Appl. Cryst. 46: 365-371.
  http://dx.doi.org/10.1107/S0021889813001295

<br /><strong>What to do in case of unsuccessful fits:</strong>
If convergence is not reached no output is generated, and only this log is 
stored in a file. On success, the resulting size distribution and data 
fit are stored to files with uncertainties.

If convergence is not reached, the following can be attempted:

1) Adjust the convergence criterion to a larger value. 
   The convergence criterion can be adjusted by the user, to support data 
   whose uncertainties are too large or too small. 
2) Verify that the parameter range of the model is appropriate. Very wide 
   ranges may prevent the correct solution from appearing within the limited
   number of iterations attempted by the program.
3) Start with a simple model. The unidirectional degeneracy of scattering 
   behaviour means that most scattering patterns can be fit using spheres. 
   Start with spheres, and move to more complex shapes afterwards. 
   Furthermore, when choosing complex models, ensure that there is evidence
   for the necessity of these complex models from alternative techniques. 
   For example, only choose "rods" when there is evidence for rods in your 
   sample (e.g. from TEM images). If a model fits your scattering pattern, 
   it does not mean it is *the* model. 

[ For more information, please see http://www.mcsas.net ]
""")

CHANGESTEXT = (u"""

Changes in v.1.2: 
- Smearing added for pinhole systems as well
- Trapezoidal and Gaussian beam profiles can be used
- Datapoints can be clipped from the beginning of the dataset
bugfixes

Changes in v.1.1:
- Slit-smearing added for Kratky-cameras
- Q-limits added

Changes in v1.0.1:
- Updated SLD range limits for models: lmadensesphere and ellipsoids
- Updated information on related literature

Changes in v1.0:
- Compiled versions available for Linux, Mac OS X and Windows.
- Histogram ranges automatically follow parameter ranges (can be disabled)
- Shannon channel estimate shown in the file dialog
- Tooltips shown when hovering over an input window
- All output stored in directories
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
        self.dataWidget = DataWidget(self)
        self.dataWidget.sigConfig.connect(self.fileWidget.setDataConfig)
        self.fileWidget.sigSelectedData.connect(self.dataWidget.onDataSelected)
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
        for settingsWidget in self.optimWidget, self.modelWidget:
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
        for settingsWidget in self.optimWidget, self.modelWidget:
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
