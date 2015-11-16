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
            return loaddatafile(fn).getDataObj()
        newIndex = len(self) # index of the next data set added
        DataList.loadData(self, sourceList = fileList, showProgress = False,
                          processSourceFunc = loaddataobj)
        self.triggerAccumulation(newIndex)

    def triggerAccumulation(self, newIndex):
        """Starts accumulation of related data sets among the currently loaded
        ones. Finally, removes such source data sets and adds the new combined
        one. newIndex is the index of the first newly loaded entries."""
        # allow accumulation of items based on the last item loaded
        lastIndex, lastData = self.currentSelection()
        if lastData is None:
            return
        avg = lastData.accumulate(self.data())
        if avg is None:
            return
        # remove the single data sets which where just loaded
        for i in range(newIndex, lastIndex + 1):
            self.setCurrentIndex(i)
            self.removeSelected()
        # add the new combined data set (again)
#        self.add(avg)
#        return
        # add the combined dls data split up per angle
        [self.add(d) for d in avg.splitPerAngle()]
        # TODO: alternatively, add a new hierarchical DataSet object containing
        # the source files as children which aren't processed further
        # FileList needs to be reconfigured towards a tree structure; while
        # being just a construtor switch the general appereance and behaviour
        # might change

    def itemDoubleClicked(self, item, column):
        valueRange = item.data().sphericalSizeEst()
        self.sigSphericalSizeRange.emit(min(valueRange), max(valueRange))
        # self.sigShannonChannels.emit(item.data().shannonChannelEst()

    def setupUi(self):
        self.listWidget.setAlternatingRowColors(True)
        setBackgroundStyleSheet(self, "./resources/background_files.svg")

    def setDataConfig(self, dataConfig):
        if self.isEmpty():
            return
        def setConfigToData(data, config = None):
            """Helper to call the appropriate method in the class hierarchy of
            the dataset."""
            data.setConfig(config)
        self.updateData(updateFunc = setConfigToData, config = dataConfig,
                        showProgress = False)
        # data = self.data(0)[0]
        # print >>sys.__stderr__, "FileList.data", type(data), id(data.config.smearing)

# vim: set ts=4 sts=4 sw=4 tw=0:
