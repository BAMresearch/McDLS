# -*- coding: utf-8 -*-
# gui/filelist.py

from __future__ import absolute_import # PEP328
from builtins import range
import logging

import os.path
from collections import OrderedDict
from gui.utils.signal import Signal
from gui.bases.datalist import DataList
from gui.utils.filedialog import getOpenFiles, getSaveDirectory
from gui.utils.translate import tr
from gui.utils import processEventLoop
from datafile import getFileFilter
from utils.lastpath import LastPath
from utils.units import ScatteringVector, ScatteringIntensity
from utils import isList
from datafile import loaddatafile
from dataobj import DataObj
from dataobj import DLSData
from getUncertainties import getUncertainties, simulate

# required for svg graphics support
from gui.liststyle import setBackgroundStyleSheet                              

from gui.qt import QtCore
from QtCore import QTime

def delay(msecs):
    dieTime = QTime.currentTime().addMSecs(msecs)
    while QTime.currentTime() < dieTime:
        processEventLoop()

class FileList(DataList):
    sigSphericalSizeRange = Signal((float, float))
    # sigShannonChannels = Signal(int) # sets warning level of "nbins"-field

    def __init__(self, *args, **kwargs):
        super(FileList, self).__init__(*args, **kwargs)
#        self.addMenuEntry(
#            name = "simulate", text = tr("add simulated"), menuStates = "*",
#            toolTip = tr("Add simulated data."),
#            callbacks = self.simulateData)

    def simulateData(self):
        print("simulateData")
        path = getSaveDirectory(self,
                    tr("Select a result directory."), LastPath.get())
        delay(300) # wait for the file dialog to go away
        if not os.path.isdir(path): return
        print("paths:", path)
        LastPath.set(os.path.dirname(path))
        ucByAngle = getUncertainties(path)
        if ucByAngle is None:
            return # nothing to do
        print("angles found:", list(ucByAngle.keys()))
        dlsData = simulate(*ucByAngle[90.])
        logging.info("Generated DLS data:\n" + str(dlsData))
        uc = dlsData.correlation.siDataU
        logging.info("Simulated DLS data with uncertainty min: {}, max: {}."
                     .format(uc.min(), uc.max()))
        DataList.loadData(self, sourceList = [dlsData], showProgress = False,
                          processSourceFunc = lambda x: x)
        print("simulateData done")

    def loadData(self, fileList = None):
        if fileList is None or type(fileList) is bool:
            fileList = getOpenFiles(self,
                # show same unit as in SASData.__init__()
                u"Load one or more .ASC data files",
                LastPath.get(), multiple = True,
                filefilter = getFileFilter()
            )
        def loaddataobj(fn):
            datafile = loaddatafile(fn)
            if datafile is None:
                return None
            dataobj = datafile.getDataObj()
            return dataobj

        config = self.configFromLast()
        nextIdx = len(self) # index of the next data set added
        # Populate the data list widget with items based on loaddataobj()
        DataList.loadData(self, sourceList = fileList, showProgress = False,
                          processSourceFunc = loaddataobj)
        # set the config of already loaded data, if any
        self.setDataConfig(config)
        # select the newly loaded data which triggers construction
        # of a DataWidget which in turn restores the last DataConfig settings
        # (if there was no data loaded yet)
        self.selectionChanged()
        self.preProcess(nextIdx)
        if config is None:
            config = self.configFromLast()
        # put the config of the last to all recently loaded
        self.setDataConfig(config)
        # update minimum uncertainty once, called too often in setConfig()
        self.updateData(updateFunc = lambda data: data.config.updateFuMin(),
                        showProgress = False)
        self.selectionChanged() # emit the last data object finally

    def configFromLast(self):
        """Get the data config of the last item in the list."""
        # FIXME for all samples and data types currently loaded,
        #       see setDataConfig()
        if self.isEmpty():
            return None
        return self.data(len(self)-1)[0].config

    def preProcess(self, firstIdx):
        """Replaces initially loaded data objects with the pre-processing
        result. firstIdx is the index of the first newly loaded entries."""
        # allow accumulation of items based on the last item loaded
        newData = self.data()[firstIdx:]
        newData = DLSData.preProcess(newData)
        if not isList(newData) or not len(newData):
            return
        # remove the single data sets which where just loaded
        self.removeItems(list(range(firstIdx, len(self))))
        # base class loadData() communicates new data via signal finally
        DataList.loadData(self, sourceList = newData,
                          showProgress = False,
                          processSourceFunc = lambda x: x)

    def itemDoubleClicked(self, item, column):
        if not hasattr(item.data(), "sphericalSizeEst"):
            return
        valueRange = item.data().sphericalSizeEst()
        self.sigSphericalSizeRange.emit(min(valueRange), max(valueRange))
        # self.sigShannonChannels.emit(item.data().shannonChannelEst()

    def setupUi(self):
        self.listWidget.setAlternatingRowColors(True)
        setBackgroundStyleSheet(self, "./resources/background_files.svg")

    def setDataConfig(self, dataConfig):
        """Propagates the given DataConfig to all DataObj in the list.
        Makes sure that all data sets have the same configuration finally.
        Disable this in order to have individual per-data-set configuration.
        """
        if self.isEmpty() or dataConfig is None:
            return
        def setConfigToData(data, config = None):
            """Helper to call the appropriate method in the class hierarchy of
            the dataset."""
            data.setConfig(config)
        self.updateData(updateFunc = setConfigToData, config = dataConfig,
                        showProgress = False)
        self.updateItems()

# vim: set ts=4 sts=4 sw=4 tw=0:
