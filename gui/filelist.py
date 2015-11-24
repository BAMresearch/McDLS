# -*- coding: utf-8 -*-
# gui/filelist.py

from __future__ import absolute_import # PEP328
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
            if datafile is not None:
                return datafile.getDataObj()
        nextIdx = len(self) # index of the next data set added
        DataList.loadData(self, sourceList = fileList, showProgress = False,
                          processSourceFunc = loaddataobj)
        self.preProcess(nextIdx)

    def preProcess(self, firstIdx):
        """Starts accumulation of related data sets among the currently loaded
        ones. Finally, removes such source data sets and adds the new combined
        one. firstIdx is the index of the first newly loaded entries."""
        # allow accumulation of items based on the last item loaded
        newData = self.data()[firstIdx:]
        if not isList(newData) or not len(newData):
            return # nothing to do
        samples = OrderedDict()
        for d in newData: # group data objects by their title
            if d.title not in samples:
                samples[d.title] = []
            samples[d.title].append(d)
        newData = [] # accumulate data objects if possible
        for title, lst in samples.iteritems():
            avg = lst[0].accumulate(lst)
            if avg is None:
                newData.extend(lst)
            else:
                newData.append(avg)
        # remove the single data sets which where just loaded
        self.removeItems(range(firstIdx, len(self)))
        # add the new combined data set (again)
#        [self.add(d) for d in newData]
#        return
        # add the combined dls data split up per angle
        for d in newData:
            try:
                [self.add(s) for s in d.splitPerAngle()]
            except AttributeError:
                self.add(d)

    def itemDoubleClicked(self, item, column):
        valueRange = item.data().sphericalSizeEst()
        self.sigSphericalSizeRange.emit(min(valueRange), max(valueRange))
        # self.sigShannonChannels.emit(item.data().shannonChannelEst()

    def setupUi(self):
        self.listWidget.setAlternatingRowColors(True)
        setBackgroundStyleSheet(self, "./resources/background_files.svg")

    def setDataConfig(self, dataConfig):
        """Propagates the given DataConfig to all DataObj in the list."""
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
