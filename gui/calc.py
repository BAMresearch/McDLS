# -*- coding: utf-8 -*-
# calc.py

import sys
import logging
import time
import os.path
import ConfigParser
import numpy
import numpy as np
from cutesnake.qt import QtCore
from QtCore import QUrl
from cutesnake.dataset import DataSet, ResultMixin
from cutesnake.utils import isList
from cutesnake.utils.lastpath import LastPath
from cutesnake.datafile import PDHFile, AsciiFile
from cutesnake.utilsgui.displayexception import DisplayException
from cutesnake.log import timestamp, addHandler
import cutesnake.log as log
from McSAS import McSAS

class SASData(DataSet, ResultMixin):
    # static settings, move this to global app settings later
    indent = "    "
    nolog = False
    # handles of gui components
    logWidget = None
    settings = None
    # McSAS algorithm instance
    mcsas = McSAS.factory()()

    @classmethod
    def load(cls, filename):
        if not os.path.isfile(filename):
            logging.warning("File '{0}' does not exist!".format(filename))
            return
        logging.info("Loading '{0}' ...".format(filename))

        sasFile = PDHFile(filename)
        sasData = cls(sasFile.name, sasFile.data)
        return sasData

    def __init__(self, title, data):
        DataSet.__init__(self, title, data)

        if self.logWidget:
            self.logWidget.clear()

        # start log file writing
        fn = self.getResultFilename("log", "this log")
        logFile = logging.FileHandler(fn, encoding = "utf8")
        oldHandler = log.log.handlers[0]
        log.addHandler(logFile)

        q, I = data[:, 0], data[:, 1]
        emin, E = 0.01, None #minimum possible error (1%)
        E = emin * I
        if data.shape[1] < 3:
            logging.warning("No error column provided! Using {}% of intensity."
                            .format(emin * 100.))
        else:
            logging.warning("Minimum uncertainty ({}% of intensity) set "
                            "for {} datapoints.".format(
                            emin * 100., sum(E > data[:, 2])))
            E = np.maximum(E, data[:, 2])

        bounds = np.array([np.pi/np.max(q), np.pi/np.min(q)])
        nreps, ncontrib, maxiter, convcrit, histbins = 3, 200, 1e4, 5., 50
        if self.settings:
            convcrit = self.settings.get("convergenceCriterion")
            nreps = self.settings.get("numReps")
            ncontrib = self.settings.get("numContribs")
            histbins = self.settings.get("histogramBins")
        logging.info("The following parameters are used for 'Analyze_1D':")
        logging.info("bounds: [{0:.4f}; {1:.4f}], Histbins: {2}"
                     .format(bounds[0], bounds[1], histbins))
        logging.info("Nsph: {0}, Nreps: {1}, MaxIter: {2}, convcrit: {3}"
                     .format(ncontrib, nreps, maxiter, convcrit))
        mcargs = dict(Emin = emin, numContribs = ncontrib, numReps = nreps,
                      histogramBins = histbins,
                      contribParamBounds = bounds,
                      maxIterations = maxiter,
                      convergenceCriterion = convcrit, doPlot = True)
        if self.nolog:
            log.removeHandler(oldHandler)
        self.mcsas.calc(Q = q, I = I, IError = E, **mcargs)
        if self.nolog:
            log.addHandler(oldHandler)

        if isList(self.mcsas.result) and len(self.mcsas.result):
            res = self.mcsas.result[0]
            if res is not None:
                self.writeDistrib(res)
                self.writeFit(res)
        else:
            logging.info("No results available!")

        self.writeSettings(mcargs)
        log.removeHandler(logFile)

    def writeFit(self, mcResult):
        self.writeResultHelper(mcResult, "fit", "fit data",
            ('fitQ', 'fitIntensityMean', 'fitIntensityStd')
        )

    def writeDistrib(self, mcResult):
        self.writeResultHelper(mcResult, "dist", "distributions",
            ('histogramXMean', 'volumeHistogramYMean',
             'volumeHistogramYStd', 'volumeHistogramMinimumRequired',
             'numberHistogramYMean', 'numberHistogramYStd',
             'numberHistogramMinimumRequired')
        )

    def writeSettings(self, mcargs):
        fn = self.getResultFilename("settings", "algorithm settings")
        config = ConfigParser.RawConfigParser()
        sectionName = "Settings"
        config.add_section(sectionName)
        for key, value in mcargs.iteritems():
            config.set(sectionName, key, value)
        config.set(sectionName, "model", self.mcsas.model.name())
        sectionName = "Model Settings"
        config.add_section(sectionName)
        for p in self.mcsas.model.params():
            config.set(sectionName, p.name(), p.value())
            config.set(sectionName, p.name()+"_min", p.min())
            config.set(sectionName, p.name()+"_max", p.max())
        with open(fn, 'w') as configfile:
            config.write(configfile)

    def getResultFilename(self, fileKey, descr):
        fn = self.getFilename(fileKey)
        logging.info("Writing {0} to:".format(descr))
        logging.info("{0}'{1}'".format(self.indent,
                                       QUrl.fromLocalFile(fn).toEncoded()))
        return fn

    def writeResultHelper(self, mcResult, fileKey, descr, columnNames):
        fn = self.getResultFilename(fileKey, descr)
        logging.info("Containing the following columns:")
        for cn in columnNames:
            logging.info("{0}[ {1} ]".format(self.indent, cn))
        data = np.vstack([mcResult[cn] for cn in columnNames]).T
        AsciiFile.writeFile(fn, data)

    def getFilename(self, kind):
        """Creates a file name from data base name, its directory and the
        current timestamp. It's created once so that all output files have
        the same base name and timestamp."""
        if not hasattr(self, "basefn") or self.basefn is None:
            self.basefn = "{0}_{1}".format(
                    os.path.join(LastPath.get(), self.title), log.timestamp())
        return "{0}_{1}.txt".format(self.basefn, kind)

def calc(filenames):
    if not isList(filenames):
        return
    if len(filenames) <= 0:
        logging.info("Please load an input file!")
        return
    for fn in filenames:
        try:
            SASData.load(fn)
        except StandardError, e:
            # on error, skip the current file
            logging.error(str(e).replace("\n"," ") + " ... skipping")
            continue
    if sys.exc_info()[0] is not None: # reraise the last error if any
        raise

# vim: set ts=4 sts=4 sw=4 tw=0:
