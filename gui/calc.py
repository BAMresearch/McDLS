# -*- coding: utf-8 -*-
# calc.py

import sys
import logging
import time
import os.path
import ConfigParser
import numpy
import numpy as np
import pickle
from cutesnake.qt import QtCore
from QtCore import QUrl
from cutesnake.dataset import DataSet, DisplayMixin
from cutesnake.utils import isList, isString, testfor
from cutesnake.utils.lastpath import LastPath
from cutesnake.datafile import PDHFile, AsciiFile
from cutesnake.utilsgui.displayexception import DisplayException
from cutesnake.log import timestamp, addHandler
import cutesnake.log as log
from mcsas.mcsas import McSAS
from utils.parameter import Histogram, RangeStats

class Calculator(object):
    _algo = None # McSAS algorithm instance
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

    def _setBaseFilename(self, dataset):
        self.basefn = "{fn}_{ts}".format(
                fn = os.path.join(LastPath.get(), dataset.title),
                ts = log.timestamp())

    def __call__(self, dataset):
        if self.model is None:
            logging.warning("No model set!")
            return
        # start log file writing
        testfor(isinstance(dataset, DataSet), StandardError,
                "{cls} requires a DataSet!".format(cls = type(self)))
        self._setBaseFilename(dataset)
        fn = self._getResultFilename("log", "this log")
        logFile = logging.FileHandler(fn, encoding = "utf8")
        oldHandler = log.log.handlers[0]
        log.addHandler(logFile)

        bounds = dataset.sphericalSizeEst()
        logging.info("The following parameters are used for 'Analyze_1D':")
        logging.info("bounds: [{0:.4f}; {1:.4f}]"
                     .format(bounds[0], bounds[1]))
        mcargs = dict(contribParamBounds = bounds,
                      doPlot = False)
        self._writeSettings(mcargs, dataset)
        if self.nolog:
            log.removeHandler(oldHandler)
        self._algo.figureTitle = os.path.basename(self.basefn)
        self._algo.calc(Q = dataset.q, I = dataset.i,
                        IError = dataset.uncertainty, **mcargs)
        if self.nolog:
            log.addHandler(oldHandler)

        if isList(self._algo.result) and len(self._algo.result):
            res = self._algo.result[0]
            # quick hack for now, will get fixed with Parameter design
            for i, p in enumerate(self.model.activeParams()):
                res = self._algo.result[i]
                self._writeDistrib(res, p.name())
                self._writeStatistics(i)
            #plotting last so calcStats is already executed.
            if res is not None:
                self._writeFit(res)
                self._writeContribs(res)
                self._algo.plot()
        else:
            logging.info("No results available!")

        log.removeHandler(logFile)

    def _writeStatistics(self, paramIndex):
        stats = dict()
        columnNames = (("lower", "upper", "weighting")
                        + RangeStats.fieldNames())
        for cn in columnNames:
            stats[cn] = []
        param = self.modelActiveParams()[paramIndex]
        param.histogram().addRange(0, numpy.inf)
        param.histogram().calcStats(paramIndex, self._algo)
        for valueRange, weighting, rangeStats in param.histogram().iterStats():
            values = valueRange + (weighting,) + rangeStats.fields
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

    def _writeDistrib(self, mcResult, paramName):
        self._writeResultHelper(mcResult,
            "dist_"+paramName, "distributions",
            ('histogramXMean', 'histogramXWidth',
             'volumeHistogramYMean', 'volumeHistogramYStd',
             'volumeHistogramMinimumRequired',
             'numberHistogramYMean', 'numberHistogramYStd',
             'numberHistogramMinimumRequired',
             'volumeHistogramYCumMean', 'volumeHistogramYCumStd',
             'numberHistogramYCumMean', 'numberHistogramYCumStd',
             ),
            extension = '.csv'
        )

    def _writeContribs(self, mcResult):
        # Writes the contribution parameters to a pickled file.
        # Can be used to continue or reanalyse a previously fitted file
        fn = self._getResultFilename("contributions",
                                     "Model contribution parameters",
                                     extension = '.pickle')
        with open(fn,'w') as fh:
            pickle.dump(mcResult['contribs'], fh)

    def _writeSettings(self, mcargs, dataset):
        if self.model is None:
            return []
        fn = self._getResultFilename("settings", "algorithm settings", 
                                     extension = '.cfg')
        config = ConfigParser.RawConfigParser()

        sectionName = "I/O Settings"
        config.add_section(sectionName)
        # the absolute filename with extension, see SASData.load()
        config.set(sectionName, 'fileName', dataset.filename)
        # the filename with timestamp of results
        config.set(sectionName, 'outputBaseName', os.path.basename(self.basefn))

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
            if p.isActive():
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
        fn = self._getResultFilename(fileKey, descr, extension = extension)
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

    def _getResultFilename(self, fileKey, descr, extension = '.txt'):
        fn = self._getFilename(fileKey, extension = extension)
        logging.info("Writing {0} to:".format(descr))
        logging.info("{0}'{1}'".format(self.indent,
                                       QUrl.fromLocalFile(fn).toEncoded()))
        return fn

    def _getFilename(self, kind, extension = '.txt'):
        """Creates a file name from data base name, its directory and the
        current timestamp. It's created once so that all output files have
        the same base name and timestamp."""
        return "{fn}_{kind}{ext}".format(
                fn = self.basefn, kind = kind, ext = extension)

# vim: set ts=4 sts=4 sw=4 tw=0:
