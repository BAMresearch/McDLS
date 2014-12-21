# -*- coding: utf-8 -*-
# gui/calc.py

from __future__ import absolute_import # PEP328
import sys
import logging
import time
import os.path
import ConfigParser
import numpy
import numpy as np
import pickle

from gui.qt import QtCore
from QtCore import QUrl
from bases.dataset import DataSet, DisplayMixin
from utils import isList, isString, testfor
from utils.lastpath import LastPath
from bases.datafile import PDHFile, AsciiFile
from gui.utils.displayexception import DisplayException
import log
from mcsas.mcsas import McSAS
from utils.parameter import Histogram, Moments, isActiveParam

class OutputFilename(object):
    """Generates output filenames with a common time stamp and logs
    appropriate messages."""
    _outDir = None   # output directory
    _basename = None # base file name for all output files
    _indent = "    "

    @property
    def basename(self):
        return self._basename

    def __init__(self, dataset):
        self._outDir = LastPath.get()
        if not os.path.isdir(self._outDir):
            logging.warning("Provided output path '{}' does not exist!"
                            .format(self._outDir))
            self._outDir = ""
        self._basename = "{title}_{ts}".format(
                title = dataset.title, ts = log.timestamp())
        # create a directory for all output files
        newDir = os.path.join(self._outDir, self._basename)
        try:
            os.mkdir(newDir)
            self._outDir = newDir
        except OSError:
            pass # on failure: no subdirectory

    def filename(self, kind = None, extension = '.txt'):
        """Creates a file name from data base name, its directory and the
        current timestamp. It's created once so that all output files have
        the same base name and timestamp."""
        fn = [self._basename]
        if isString(kind):
            fn += ["_", kind]
        if isString(extension):
            fn += extension
        return os.path.join(self._outDir, "".join(fn))

    def filenameVerbose(self, kind, descr, extension = '.txt'):
        """Returns the file name as in filename() and logs a descriptive
        message containing the full file name which is usually click-able."""
        fn = self.filename(kind, extension = extension)
        logging.info("Writing {0} to:".format(descr))
        logging.info("{0}'{1}'".format(self._indent,
                                       QUrl.fromLocalFile(fn).toEncoded()))
        return fn

class Calculator(object):
    _algo = None  # McSAS algorithm instance
    _outFn = None # handles output file names, creates directories
    # static settings, move this to global app settings later
    indent = "    "
    nolog = False

    def __init__(self):
        self._algo = McSAS.factory()()

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

    def __call__(self, dataset):
        if self.model is None:
            logging.warning("No model set!")
            return
        # start log file writing
        testfor(isinstance(dataset, DataSet), StandardError,
                "{cls} requires a DataSet!".format(cls = type(self)))
        self._outFn = OutputFilename(dataset)
        fn = self._outFn.filenameVerbose("log", "this log")
        logFile = logging.FileHandler(fn, encoding = "utf8")
        oldHandler = log.log.handlers[0]
        log.addHandler(logFile)

        # obsolete? 2014-11-28:
        # bounds = dataset.sphericalSizeEst()
        # logging.info("The following parameters are used for 'Analyze_1D':")
        # logging.info("bounds: [{0:.4f}; {1:.4f}]"
        #              .format(bounds[0], bounds[1]))
        mcargs = dict(doPlot = False)
        self._writeSettings(mcargs, dataset)
        if self.nolog:
            log.removeHandler(oldHandler)
        self._algo.calc(SASData = dataset, **mcargs)
        if self.nolog:
            log.addHandler(oldHandler)

        if isList(self._algo.result) and len(self._algo.result):
            res = self._algo.result[0]
            # quick hack for now, will get fixed with Parameter design
            for i, p in enumerate(self.model.activeParams()):
                self._writeDistrib(p)
                self._writeStatistics(i, p)
            # plotting last so calcStats is already executed.
            if res is not None:
                self._writeFit(res)
                self._writeContribs(res)
                self._algo.plot(outputFilename = self._outFn)
        else:
            logging.info("No results available!")

        log.removeHandler(logFile)

    def _writeStatistics(self, paramIndex, param):
        stats = dict()
        columnNames = (("lower", "upper", "weighting")
                        + Moments.fieldNames())
        for cn in columnNames:
            stats[cn] = []
        # not optimal, but for now, it helps
        for h in param.histograms():
            values = h.xrange + (h.yweight,) + h.moments.fields
            for name, value in zip(columnNames, values):
                if stats.get(name) is None:
                    stats[name] = []
                stats[name].append(AsciiFile.formatValue(value))
        self._writeResultHelper(stats, "stats_"+param.name(),
                                "distribution statistics",
                                columnNames, extension = '.csv')

    def _writeFit(self, mcResult):
        self._writeResultHelper(mcResult, "fit", "fit data",
            ('fitQ', 'fitIntensityMean', 'fitIntensityStd'),
            extension = '.csv'
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
                extension = '.csv'
            )

    def _writeContribs(self, mcResult):
        # Writes the contribution parameters to a pickled file.
        # Can be used to continue or reanalyse a previously fitted file
        fn = self._outFn.filenameVerbose("contributions",
                                         "Model contribution parameters",
                                         extension = '.pickle')
        with open(fn,'w') as fh:
            pickle.dump(mcResult['contribs'], fh)

    def _writeSettings(self, mcargs, dataset):
        if self.model is None:
            return []
        fn = self._outFn.filenameVerbose("settings", "algorithm settings",
                                         extension = '.cfg')
        config = ConfigParser.RawConfigParser()

        sectionName = "I/O Settings"
        config.add_section(sectionName)
        # the absolute filename with extension, see SASData.load()
        config.set(sectionName, 'fileName', dataset.filename)
        # the filename with timestamp of results
        config.set(sectionName, 'outputBaseName', self._outFn.basename)

        sectionName = "MCSAS Settings"
        config.add_section(sectionName)
        for key, value in mcargs.iteritems():
            config.set(sectionName, key, value)
        for p in self.algo.params():
            config.set(sectionName, p.name(), p.value())
        config.set(sectionName, "model", self.model.name())
        # We don't have to do anything with these yet, but storing them for now:
        config.set(sectionName, "Q limits", 
                np.array([np.min(dataset.q),np.max(dataset.q)]) )

        sectionName = "Model Settings"
        config.add_section(sectionName)
        for p in self.modelParams():
            if isActiveParam(p):
                config.set(sectionName, p.name()+"_min", p.min())
                config.set(sectionName, p.name()+"_max", p.max())
            else:
                config.set(sectionName, p.name(), p.value())
        with open(fn, 'w') as configfile:
            config.write(configfile)

    def _writeResultHelper(self, mcResult, fileKey, descr, columnNames,
                           extension = '.txt'):
        if not all(cn in mcResult for cn in columnNames):
            logging.warning('Result does not contain the requested data')
            return
        fn = self._outFn.filenameVerbose(fileKey, descr, extension = extension)
        logging.info("Containing the following columns:")
        cwidth = max([len(cn) for cn in columnNames])
        fmt = "{0}[ {1:" + str(cwidth) + "s} ]"
        for cn in columnNames:
            msg = fmt.format(self.indent, cn)
            peek = np.ravel(mcResult[cn])
            if len(peek) < 3:
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
