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
            if datafile is None:
                return None
            dataobj = datafile.getDataObj()
            return dataobj

        nextIdx = len(self) # index of the next data set added
        DataList.loadData(self, sourceList = fileList, showProgress = False,
                          processSourceFunc = loaddataobj)
        self.postProcess(nextIdx)
        # put the config of the last to all recently loaded
        self.setDataConfig(self.configFromLast())
        self.configFromLast().overrideDefaults()
        self.selectionChanged()

    def configFromLast(self):
        """Get the data config of the last item in the list."""
        if self.isEmpty():
            return None
        return self.data(len(self)-1)[0].config

    def postProcess(self, firstIdx):
        """Starts accumulation of related data sets among the currently loaded
        ones. Finally, removes such source data sets and adds the new combined
        one. firstIdx is the index of the first newly loaded entries."""
        # allow accumulation of items based on the last item loaded
        newData = self.data()[firstIdx:]
        if not isList(newData) or not len(newData):
            return # nothing to do
        samples = OrderedDict()
        def makeKey(data):
            key = (data.title,)
            if hasattr(data, "angles"):
                key += tuple(data.angles)
            return key
        for d in newData: # group data objects by their title
            key = makeKey(d)
            if key not in samples:
                samples[key] = []
            samples[key].append(d)
        newData = [] # accumulate data objects if possible
        for dummy, lst in samples.iteritems():
            avg = lst[0].accumulate(lst)
            if avg is None:
                newData.extend(lst)
            else:
                newData.append(avg)
        # remove the single data sets which where just loaded
        self.removeItems(range(firstIdx, len(self)))
        # add the combined dls data split up per angle
        def splitUp(d):
            try: # perhaps test for isinstance(d, DLSData) instead
                return d.splitPerAngle()
            except AttributeError:
                return (d,)
        newData = [s for dl in (splitUp(d) for d in newData) for s in dl]
        # use loadData() to communicate new data finally via signal
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
            if data.sampleName != config.sampleName:
                return
            data.setConfig(config)
        self.updateData(updateFunc = setConfigToData, config = dataConfig,
                        showProgress = False)
        self.selectionChanged()

# vim: set ts=4 sts=4 sw=4 tw=0:
