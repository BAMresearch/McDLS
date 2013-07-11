# -*- coding: utf-8 -*-
# calc.py

import logging
import os.path
from cutesnake.dataset import DataSet, ResultMixin
from cutesnake.utils import isList
from cutesnake.utils.lastpath import LastPath
from cutesnake.datafile import PDHFile, AsciiFile
from cutesnake.utilsgui.displayexception import DisplayException
from cutesnake.log import Log

import mcsas
from mcsas.McSAS import ScatteringModel

import numpy as np

import time
from cutesnake.utilsgui import processEventLoop

import numpy

class SASData(DataSet, ResultMixin):
    logWidget = None
    settings = None
    mcsas = mcsas.McSAS.factory()()

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

        q, I = data[:, 0], data[:, 1]
        emin, E = 0.01, None
        if data.shape[1] < 3:
            factor = 0.01
            E = factor * I
            logging.warning("No error column provided! Using {0}% of intensity."
                            .format(factor * 100.))
        else:
            E = np.maximum(emin * I, data[:, 2])

        bounds = np.array([np.pi/np.max(q), np.pi/np.min(q)])
        nreps, ncontrib, maxiter, convcrit, histbins = 3, 200, 1e4, 1.3, 50
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
        if (not isinstance(self.mcsas.model, ScatteringModel) and
                issubclass(self.mcsas.model, ScatteringModel)):
            self.mcsas.model = self.mcsas.model()
        self.mcsas.calc(Q = q, I = I, IError = E, **mcargs)

        if len(self.mcsas.result):
            res = self.mcsas.result[0]
            if res is not None:
                self.writeDistrib(res, mcargs)
                self.writeFit(res)
        else:
            logging.info("No results available!")
        self.writeLog()

    def writeFit(self, mcResult):
        fn = self.getFilename("fit")
        logging.info("Writing fit to '{0}'".format(fn))
        data = np.vstack((mcResult['fitQ'],
                          mcResult['fitIntensityMean'],
                          mcResult['fitIntensityStd'])).T
        AsciiFile.writeFile(fn, data)

    def writeDistrib(self, mcResult, mcargs):
        fn = self.getFilename("dist")
        logging.info("Writing distribution to '{0}'".format(fn))
        data = np.vstack((mcResult['histogramXMean'],
                          mcResult['volumeHistogramYMean'],
                          mcResult['volumeHistogramYStd'],
                          mcResult['volumeHistogramMinimumRequired'],
                          mcResult['numberHistogramYMean'],
                          mcResult['numberHistogramYStd'],
                          mcResult['numberHistogramMinimumRequired'])).T
        AsciiFile.writeFile(fn, data)
        import ConfigParser
        config = ConfigParser.RawConfigParser()
        fn = self.getFilename("setting")
        sectionName = "Settings"
        config.add_section(sectionName)
        for key, value in mcargs.iteritems():
            config.set(sectionName, key, value)
        config.set(sectionName, "model", self.mcsas.model.name)
        sectionName = "Model Settings"
        config.add_section(sectionName)
        for ptype in self.mcsas.model.parameters:
            p = getattr(self.mcsas.model, ptype.name())
            config.set(sectionName, p.name(), p.value())
            config.set(sectionName, p.name()+"_min", p.min())
            config.set(sectionName, p.name()+"_max", p.max())
        with open(fn, 'wb') as configfile:
            config.write(configfile)

    def writeLog(self):
        if not self.logWidget:
            return
        fn = self.getFilename("log")
        logging.info("Writing log to '{0}'".format(fn))
        self.logWidget.saveToFile(filename = fn)

    def getFilename(self, kind):
        if not hasattr(self, "basefn") or self.basefn is None:
            self.basefn = "{0}_{1}".format(
                    os.path.join(LastPath.get(), self.title), Log.timestamp())
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
            DisplayException(e)
            logging.error(str(e).replace("\n"," ") + " ... skipping")
            continue

# vim: set ts=4 sts=4 sw=4 tw=0:
