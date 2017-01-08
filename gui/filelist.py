# -*- coding: utf-8 -*-
# gui/filelist.py

from __future__ import absolute_import # PEP328
from builtins import range
import logging

from collections import OrderedDict
from gui.utils.signal import Signal
from gui.bases.datalist import DataList
from gui.utils.filedialog import getOpenFiles
from datafile import getFileFilter
from utils.lastpath import LastPath
from utils.units import ScatteringVector, ScatteringIntensity
from utils import isList
from datafile import loaddatafile
from dataobj import DataObj
from dataobj import DLSData

# required for svg graphics support
from gui.liststyle import setBackgroundStyleSheet                              

class FileList(DataList):
    sigSphericalSizeRange = Signal((float, float))
    # sigShannonChannels = Signal(int) # sets warning level of "nbins"-field

    def loadData(self, fileList = None):
        if fileList is None or type(fileList) is bool:
            fileList = getOpenFiles(self,
                # show same unit as in SASData.__init__()
                u"Load one or more data files with q({qu}) and intensity({iu})"
                .format(qu = ScatteringVector(u"nm⁻¹").displayMagnitudeName,
                        iu = ScatteringIntensity(u"(m sr)⁻¹").displayMagnitudeName),
                LastPath.get(), multiple = True,
                filefilter = getFileFilter()
            )
        def loaddataobj(fn):
            datafile = loaddatafile(fn)
            if datafile is None:
                return None
            dataobj = datafile.getDataObj()
            return dataobj

        nextIdx = len(self) # index of the next data set added
        # Populate the data list widget with items based on loaddataobj()
        DataList.loadData(self, sourceList = fileList, showProgress = False,
                          processSourceFunc = loaddataobj)
        # select the newly loaded data which triggers construction
        # of a DataWidget which in turn restores the last DataConfig settings
        self.selectionChanged()
        self.preProcess(nextIdx)
        # put the config of the last to all recently loaded
        config = self.configFromLast()
        self.setDataConfig(config)
        if config is not None:
            config.overrideDefaults()
        self.selectionChanged()

    def configFromLast(self):
        """Get the data config of the last item in the list."""
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
        if self.isEmpty():
            return
        def setConfigToData(data, config = None):
            """Helper to call the appropriate method in the class hierarchy of
            the dataset."""
            # set this only for equal sample names,
            # -> None == None for SASData (no sample name yet)
            if data.sampleName != config.sampleName:
                return
            data.setConfig(config)
        self.updateData(updateFunc = setConfigToData, config = dataConfig,
                        showProgress = False)

# vim: set ts=4 sts=4 sw=4 tw=0:
