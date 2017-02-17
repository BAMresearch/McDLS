# -*- coding: utf-8 -*-
# gui/calc.py

from __future__ import absolute_import # PEP328
from future import standard_library
standard_library.install_aliases()
from builtins import zip
from builtins import str
from builtins import object
import sys
import logging
import time
import os.path
import codecs
from collections import OrderedDict
try: 
    import configparser
except ImportError: 
    import ConfigParser as configparser
import numpy as np
import pickle

from gui.qt import QtCore
from QtCore import QUrl
from bases.dataset import DataSet
from utils import (isList, isString, isNumber, testfor, isMac, fixFilename,
                   mcopen)
from utils.lastpath import LastPath
from utils.units import Angle
from datafile import PDHFile, AsciiFile
from gui.utils.displayexception import DisplayException
import log
from mcsas.mcsas import McSAS
from utils.parameter import Histogram, Moments, isActiveFitParam
from dataobj import SASData
from utils.hdf import HDFMixin
from mcsas.plotting import PlotSeriesStats

DEFAULTSECT = configparser.DEFAULTSECT

def cfgwrite(self, fp):
    """Write an .ini-format representation of the configuration state."""
    if self._defaults:
        fp.write("[%s]\n" % DEFAULTSECT)
        for (key, value) in list(self._defaults.items()):
            fp.write("%s = %s\n" % (key, str(value).replace('\n', '\n\t')))
        fp.write("\n")
    for section in self._sections:
        fp.write("[%s]\n" % section)
        for (key, value) in list(self._sections[section].items()):
            if key == "__name__":
                continue
            if (value is not None) or (self._optcre == self.OPTCRE):
                key = " = ".join((key, str(value).replace('\n', '\n\t')))
            fp.write("%s\n" % (key))
        fp.write("\n")

configparser.RawConfigParser.write = cfgwrite

class OutputFilename(object):
    """Generates output filenames with a common time stamp and logs
    appropriate messages."""
    _outDir = None   # output directory
    _basename = None # base file name for all output files
    _indent = "    "

    @property
    def basename(self):
        return self._basename

    @property
    def outDir(self):
        return self._outDir

    def __init__(self, dataset, createDir = True):
        self._outDir = LastPath.get()
        if not os.path.isdir(self._outDir):
            logging.warning("Output path '{}' does not exist!"
                            .format(self._outDir))
            self._outDir = ""
        self._basename = u"{title} {ts}".format(
                title = dataset.title, ts = log.timestamp())
        if hasattr(dataset, "anglesToStr"):
            self._basename += u" [{}]".format(dataset.anglesToStr)
        if hasattr(dataset, "measIndicesStr"):
            self._basename += u" ({})".format(dataset.measIndicesStr)
        if not createDir:
            return
        # create a directory for all output files
        newDir = os.path.join(self._outDir, self._basename)
        try:
            os.mkdir(newDir)
            self._outDir = newDir
        except OSError:
            logging.warning("Failed to create directory '{}'!".format(newDir))
            pass # on failure: no subdirectory

    def filename(self, kind = None, extension = '.txt'):
        """Creates a file name from data base name, its directory and the
        current timestamp. It's created once so that all output files have
        the same base name and timestamp."""
        fn = [self._basename]
        if isString(kind) and len(kind):
            fn += ["_", kind]
        if isString(extension):
            fn += extension
        fn = os.path.join(self._outDir, "".join(fn))
        return fixFilename(fn)

    def filenameVerbose(self, kind, descr, extension = '.txt'):
        """Returns the file name as in filename() and logs a descriptive
        message containing the full file name which is usually click-able."""
        fn = self.filename(kind, extension = extension)
        logging.info("Writing {0} to:".format(descr))
        fnUrl = fn
        if fnUrl.startswith('\\\\?'):
            fnUrl = fnUrl[4:]
        logging.info("{0}'{1}'".format(self._indent,
            QUrl.fromLocalFile(fnUrl).toEncoded()))
        return fn

class Calculator(HDFMixin):
    _algo = None  # McSAS algorithm instance
    _outFn = None # handles output file names, creates directories
    _series = None # stores results of multiple data sets for a final summary
    # static settings, move this to global app settings later
    indent = "    "
    nolog = False

    def __init__(self):
        self._algo = McSAS.factory()()

    def hdfWrite(self, hdf):
        """ write a calculator configuration. """
        hdf.writeMember(self, "algo")
        hdf.writeMember(self.algo, "data")
        hdf.writeMember(self, "model")
        # for p in self.model.params():
        #     logging.debug("Writing model parameter: {} value: {} to HDF5".format(p.name(), p.value()))
        #     hdf.writeMember(self.model, p.name())

    def hdfLoad(self, filehandle):
        """ load a calculator configuration """
        pass

    @property
    def algo(self):
        return self._algo

    @property
    def model(self):
        return self._algo.model

    @model.setter
    def model(self, newModel):
        self._algo.model = newModel

    def modelParams(self):
        if self.model is None:
            return []
        # access the property to change it permanently
        return [getattr(self.model, p.name()) for p in self.model.params()]

    def modelActiveParams(self):
        if self.model is None:
            return []
        # access the property to change it permanently
        return [getattr(self.model, p.name()) for p in self.model.activeParams()]

    def stop(self):
        self._algo.stop = True

    def isStopped(self):
        return self._algo.stop

    def prepare(self):
        """Resets series data. Supposed to be called before each run of
        multiple __call__() invokations."""
        self._series = OrderedDict()

    def __call__(self, dataset):
        if self.model is None:
            logging.warning("No model set!")
            return
        # start log file writing
        testfor(isinstance(dataset, DataSet), Exception,
                "{cls} requires a DataSet!".format(cls = type(self)))
        self._outFn = OutputFilename(dataset)
        fn = self._outFn.filenameVerbose("log", "this log")
        logFile = logging.FileHandler(fn, encoding = "utf8")
        widgetHandler = log.getWidgetHandlers()[0]
        log.replaceHandler(widgetHandler) # remove everything else
        log.addHandler(logFile)
        try:
            # show last lines about pdf output from separate thread
            widgetHandler.widget.addWatchDir(self._outFn.outDir)
        except AttributeError:
            pass

        self._writeSettings(dict(), dataset)
        if self.nolog: # refers to the widgethandler
            log.removeHandler(widgetHandler)
        #set data in the algorithm
        self._algo.data = dataset

        # write HDF5, show exceptions traceback, if any
        try:
            self.hdfStore(self._outFn.filenameVerbose(
                "hdf5archive", "Complete state of the calculation",
                extension = '.mh5'), rootLocation = "mcsasentry")
        except Exception as e:
            import traceback
            print(traceback.format_exc())

        self._algo.calc()
        if self.nolog:
            log.addHandler(widgetHandler)

        if isList(self._algo.result) and len(self._algo.result):
            res = self._algo.result[0]
            # quick hack for now, will get fixed with Parameter design
            for p in self.model.activeParams():
                self._writeDistrib(p)
                self._writeStatistics(p)
            self._updateSeries(dataset, self.model)
            # plotting last so stats were already calculated
            if res is not None:
                self._writeFit(res)
                self._writeContribs(res)
                self._algo.plot(outputFilename = self._outFn,
                                autoClose = self._algo.autoClose())
        else:
            logging.info("No results available!")

        log.removeHandler(logFile)

    def _updateSeries(self, data, model):
        if not self.algo.seriesStats():
            return
        def addSeriesData(uid, hist, seriesKey):
            if uid not in self._series:
                self._series[uid] = []
            self._series[uid].append((seriesKey, hist.moments.fields))
        def makeId(data, hist):
            """Derive a unique identifier for each pair of sample and
            histogram."""
            uid = (data.sampleName, data.seriesKeyName,
                   (hist.param.name(),)
                   + tuple(hist.param.unit().toDisplay(x) for x in h.xrange)
                   + (h.yweight,))
            return uid

        for p in model.activeParams():
            for h in p.histograms():
                addSeriesData(makeId(data, h), h, data.seriesKeyValue)

    def postProcess(self):
        if not self.algo.seriesStats():
            return
        # works similar to _writeStatistics() but not using parameters

        class DummyDataSet(object):
            """Just for the file name formatting."""
            title = u"series statistics"

        def processSeries(series):
            seriesPlot = PlotSeriesStats()
            # data formatted for file output, gathered across histograms
            fileData = dict()
            columnNames = ( # columns appearing in the output file header
                ["seriesKey", "param", "lower", "upper", "weighting"]
                + list(Moments.fieldNames()))
            for seriesItem in series.items():
                processSeriesStats(seriesItem, seriesPlot, fileData, columnNames)
            # since we are the last writer, changing outFn doesn't hurt
            self._outFn = OutputFilename(DummyDataSet, createDir = False)
            self._writeResultHelper(fileData, "", "series statistics",
                                    columnNames, extension = '.dat')
            seriesPlot.show()

        def processSeriesStats(seriesItem, seriesPlot, fileData, columnNames):
            # gather data values indexed by columns names first
            stats = dict()
            ((sampleName, seriesKeyName, (pname, lo, hi, weight)),
              valuePairs) = seriesItem
            columnNames[0] = seriesKeyName.replace(" ", "_")
            for seriesKey, moments in valuePairs:
                values = (seriesKey, pname, lo, hi, weight) + moments
                for name, value in zip(columnNames, values):
                    if name not in stats:
                        stats[name] = []
                    # for plotting below, no float-str conversion here
                    stats[name].append(value)
            # convert numerical stats to proper formatted text for file output
            for key, values in stats.items():
                # proper float-str formatting for text file output
                if key not in fileData:
                    fileData[key] = []
                for value in values:
                    if isList(value):
                        value = ";".join([AsciiFile.formatValue(v)
                                          for v in value])
                    fileData[key].append(AsciiFile.formatValue(value))
            # simple statistics plotting, kind of a prototype for now ...
            stats["seriesKeyName"] = seriesKeyName
            stats["seriesKey"] = stats[columnNames[0]]
            stats["cfg"] = u"{param} [{lo},{hi}] {w}".format(
                            param = pname, lo = lo, hi = hi, w = weight)
            stats["title"] = sampleName
            seriesPlot.plot(stats)

        if isMac():
            processSeries(self._series)
        else:
            from multiprocessing import Process
            proc = Process(target = processSeries, args = (self._series,))
            proc.start()

    def _writeStatistics(self, param):
        """Gathers the statistics column-wise first and converts them to a
        row oriented text file in _writeResultHelper()"""
        stats = dict()
        columnNames = (("lower", "upper", "weighting")
                        + Moments.fieldNames())
        # not optimal, but for now, it helps
        for h in param.histograms():
            values = h.xrange + (h.yweight,) + h.moments.fields
            for name, value in zip(columnNames, values):
                if name not in stats:
                    stats[name] = []
                stats[name].append(AsciiFile.formatValue(value))
        self._writeResultHelper(stats, "stats_"+param.name(),
                                "distribution statistics",
                                columnNames, extension = '.dat')

    def _writeFit(self, mcResult):
        columnNames = ('fitX0', 'dataMean', 'dataStd',
                       'fitMeasValMean', 'fitMeasValStd')
        self._writeResultHelper(mcResult, "fit", "fit data",
            columnNames,
            extension = '.dat'
        )

    def _writeDistrib(self, param):
        for h in param.histograms():
            histRes = dict(
                xMean = h.xMean, xWidth = h.xWidth,
                yMean = h.bins.mean, yStd = h.bins.std,
                Obs = h.observability,
                cdfMean = h.cdf.mean, cdfStd = h.cdf.std
            )
            self._writeResultHelper(histRes, str(h), "distributions",
                ("xMean", "xWidth", # fixed order of result columns
                 "yMean", "yStd", "Obs",
                 "cdfMean", "cdfStd"),
                extension = '.dat'
            )

    def _writeContribs(self, mcResult):
        # Writes the contribution parameters to a pickled file.
        # Can be used to continue or reanalyse a previously fitted file
        fn = self._outFn.filenameVerbose("contributions",
                                         "Model contribution parameters",
                                         extension = '.pickle')
        with mcopen(fn, 'wb') as fh:
            pickle.dump(mcResult['contribs'], fh)

    def _writeSettings(self, mcargs, dataset):
        if self.model is None:
            return []
        fn = self._outFn.filenameVerbose("settings", "algorithm settings",
                                         extension = '.cfg')
        config = configparser.RawConfigParser()

        sectionName = "I/O Settings"
        config.add_section(sectionName)
        # the absolute filename with extension, see SASData.load()
        config.set(sectionName, 'fileName', dataset.filename)
        # the filename with timestamp of results
        config.set(sectionName, 'outputBaseName', self._outFn.basename)

        sectionName = "MCSAS Settings"
        config.add_section(sectionName)
        for key, value in mcargs.items():
            config.set(sectionName, key, value)
        for p in self.algo.params():
            config.set(sectionName, p.name(), p.value())
        config.set(sectionName, "model", self.model.name())
        # We don't have to do anything with these yet, but storing them for now:
        if isinstance(dataset, SASData): # useful with SAS data only
            config.set(sectionName, "X0 limits", str(dataset.x0.limit))

        sectionName = "Model Settings"
        config.add_section(sectionName)
        for p in self.modelParams():
            if isActiveFitParam(p):
                config.set(sectionName, p.name()+"_min", p.min())
                config.set(sectionName, p.name()+"_max", p.max())
            else:
                config.set(sectionName, p.name(), p.value())
        with mcopen(fn, 'w') as configfile:
            config.write(configfile)

    def _writeResultHelper(self, mcResult, fileKey, descr, columnNames,
                           extension = '.txt'):
        if not all(cn in mcResult for cn in columnNames):
            logging.warning(
                'Writing results: some of the requested {d} not found!'
                .format(d = descr))
            logging.warning(str(columnNames))
            return
        fn = self._outFn.filenameVerbose(fileKey, descr, extension = extension)
        logging.info("Containing the following columns:")
        cwidth = max([len(cn) for cn in columnNames])
        fmt = "{0}[ {1:" + str(cwidth) + "s} ]"
        for cn in columnNames:
            msg = fmt.format(self.indent, cn)
            peek = np.ravel(mcResult[cn])
            if len(peek) > 2:
                continue
            for value in peek[0:2]:
                try:
                    msg += " {0: .4e}".format(value)
                except ValueError:
                    msg += " {0: >14s}".format(value)
            logging.info(msg)
        # write header:
        AsciiFile.writeHeaderLine(fn, columnNames)
        # (additional header lines can be appended if necessary)
        data = np.vstack([mcResult[cn] for cn in columnNames]).T
        # append to the header, do not overwrite:
        AsciiFile.appendFile(fn, data)

# vim: set ts=4 sts=4 sw=4 tw=0:
