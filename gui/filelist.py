# -*- coding: utf-8 -*-
# gui/filelist.py

from __future__ import absolute_import # PEP328
import logging

from gui.utils.signal import Signal
from gui.bases.datalist import DataList
from gui.utils.filedialog import getOpenFiles
from datafile import getFileFilter
from utils.lastpath import LastPath
from utils.units import ScatteringVector, ScatteringIntensity
from datafile import loaddatafile
from dataobj import DataObj

# required for svg graphics support
from gui.liststyle import setBackgroundStyleSheet                              
import sys

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
        # populates to data list widget with items based on the return of
        # processSourceFunc(filename)
        def loaddataobj(fn):
            datafile = loaddatafile(fn)
            if datafile is None:
                return None
            dataobj = datafile.getDataObj()
            # set the config of the last item by default
            dataobj.setConfig(self.configFromLast())
            return dataobj

        DataList.loadData(self, sourceList = fileList, showProgress = False,
                          processSourceFunc = loaddataobj)

    def configFromLast(self):
        """Get the data config of the last item in the list."""
        if self.isEmpty():
            return None
        return self.data(len(self)-1)[0].config

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
            data.setConfig(config)
        self.updateData(updateFunc = setConfigToData, config = dataConfig,
                        showProgress = False)
        self.selectionChanged()

# vim: set ts=4 sts=4 sw=4 tw=0:
