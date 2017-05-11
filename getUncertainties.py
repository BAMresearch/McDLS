#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# retrieves measurement uncertainties from McDLS result files and directories
# which are possibly filtered, cropped and averaged by scattering angle

from __future__ import absolute_import, unicode_literals
from future import standard_library
standard_library.install_aliases()
from builtins import str
import sys
import os
import re
from collections import OrderedDict
import numpy as np
import logging
logging.basicConfig(level = logging.INFO)
#np.seterr(all='raise', under = 'ignore')

from utils import isString, isList
from utils.units import Deg

def angleFromFilename(name):
    match = re.findall('\[([^°\]]+)°\]', name)
    try:
        return float(match[0])
    except (IndexError, ValueError):
        return None

def findFitOutput(dataPath):
    if not os.path.isdir(dataPath):
        logging.error("Provided path '{}' is not a directory!"
                      .format(dataPath))
        return
    fitFiles = []
    for dirpath, dirnames, filenames in os.walk(dataPath):
        goodFiles = []
        for fn in filenames:
            if not fn.endswith("_fit.dat"):
                continue
            goodFiles.append((dirpath, fn, angleFromFilename(fn)))
        fitFiles += goodFiles
    return fitFiles

def processFitFile(dirpath, filename, *args):
    logging.info("processing '{}'".format(os.path.join(dirpath, filename)))
    filepath = os.path.join(dirpath, filename)
    lines = None
    with open(filepath, 'rb') as fd:
        lines = [l.decode("utf8") for l in fd.readlines()]
    if lines is None or len(lines) < 2:
        logging.error("No lines found!")
        return
    fields = lines[0].split()
    colidx = [i for i, field in enumerate(fields) if "dataStd" in field]
    if not len(colidx) or not "fitX0" in fields[0]:
        logging.error("No uncertainty data found! Fields found: {}".format(fields))
        return
    colidx = colidx[0]
    xvec, yvec = [], []
    for l in lines[1:]:
        fields = l.split()
        try:
            x, y = float(fields[0]), float(fields[colidx])
            xvec.append(x)
            yvec.append(y)
        except ValueError:
            continue
    if len(xvec) != len(yvec):
        logging.error("Could not extract all uncertainty values, "
                      "skipping this one!")
        return
    return np.array((xvec, yvec)).T

import matplotlib
import matplotlib.pyplot as pyplot
from mcsas.plotting import PlotResults

class PlotUncertainty(object):
    _figure = None
    _axes = None
    _cmap = None # colormap
    _idx = None  # number of plots, counting
    _suffix = None # title suffix

    @property
    def titleFigure(self):
        title = "Uncertainties"
        if isString(self._suffix):
            title = "{} {}".format(title, self._suffix)
        return title

    @property
    def titleAxes(self): return "extracted " + self.titleFigure.lower()

    @property
    def scalex(self): return "log"

    @property
    def labelx(self): return "tau"

    @property
    def labely(self):
        lbl = "Uncertainty"
        if isString(self._suffix):
            lbl = "{} {}".format(lbl, self._suffix)
        return lbl

    def __init__(self, count = 0, suffix = None):
        self._suffix = suffix
        self._figure = pyplot.figure(figsize = (15, 8), dpi = 80,
                              facecolor = 'w', edgecolor = 'k')
        self._figure.canvas.set_window_title(self.titleFigure)
        self._axes = pyplot.subplot(xscale = self.scalex)
        self._axes.set_title(self.titleAxes)
        self._axes.set_xlabel(self.labelx)
        self._axes.set_ylabel(self.labely)
        if count > 0:
            # https://matplotlib.org/examples/color/colormaps_reference.html
            self._cmap = pyplot.cm.get_cmap("gist_ncar", count+1)
        self._idx = 0

    def plot(self, uc, label):
        if not isinstance(uc, np.ndarray):
            return
#        logging.info("plotting {} {}".format(type(uc), uc.shape))
        kwargs = dict(label = label)
        if self._cmap is not None:
            kwargs['color'] = self._cmap(self._idx)
        self._axes.plot(uc[:,0], uc[:,1], **kwargs)
        self._idx += 1

    def show(self):
        self._axes.legend(loc = 1, fancybox = True)
        pyplot.setp(pyplot.gca().get_legend().get_texts(), fontsize = 'small')
        PlotResults.plotGrid(self._axes)
#        pyplot.show()

def launchPlot(uc):
    from multiprocessing import Process, Queue
    proc = Process(target = plot, args = [uc])
    proc.start()

def simulate(uc, doPlot = True):
    try:
        assert len(uc.shape) == 2
    except (AttributeError, AssertionError):
        return None
    from models.dlssphere import DLSSphere
    from dataobj import DLSData
    from utils.units import K, Vis, NM, Deg, Sec, MSec
    # set up the dummy DataObj
    dlsData = DLSData(title = "simulated title")
    dlsData.setSampleName("simulated sample")
    dlsData.setFilename("simulated file")
    # usual values from measurements
    dlsData.setTemperature(K.toSi(296))
    dlsData.setViscosity(Vis.toSi(0.932))
    dlsData.setRefractiveIndex(1.33200)
    dlsData.setWavelength(NM.toSi(632.8))
    # init data for a single angle
    dlsData.setAngles(Deg, Deg.toSi(np.array((90.,))))
    dlsData.setTau(MSec, MSec.toDisplay(uc[:, 0]))
    dlsData.setCorrelation(np.zeros_like(uc[:,1]).reshape((-1, 1)))
    dlsData.initConfig()
    # set up the DLS model
    model = DLSSphere()
#    model.radius.setActive(False)
#    model.radius.setDisplayValue(50.)
    # calculate finally
#    from scipy.stats import norm
#    sizes = numpy.arange(200)
#    distrib = numpy.array([norm.pdf(x, 50.) for x in sizes])
#    sizes[distrib * sizes > 0.01]
    contribs = np.array([NM.toSi(x) for x in (45., 50., 50., 55.)])
    modelData = model.calc(dlsData, contribs.reshape((-1, 1)),
                           compensationExponent = 1.)
    dlsData.setCorrelation(modelData.chisqrInt.reshape((-1, 1)), uc[:,1].reshape((-1, 1)))
    return dlsData

def dummy(doPlot):
    if doPlot:
        from multiprocessing import Process, Queue
        # multithreaded plotting also logs to file
    #    pkwargs["logToFile"] = True
    #    q = Queue() # allow communication between processes
    #    pkwargs["queue"] = q
        # set up the result data, same as in the McSAS class
        result = dict(fitX0 = dlsData.x0.binnedData,
                      fitMeasValMean = np.zeros_like(uc[:,1]).reshape((1, -1)),
                      times = None)
        proc = Process(target = PlotResults, args = ((result,), dlsData))
        proc.start()
    from datafile.cgsfile import CGSFile
    from gui.calc import OutputFilename
    outFn = OutputFilename(dlsData, createDir = False)
    fn = outFn.filenameVerbose("sim", "simulated DLS data", extension = ".ASC")
    fn = os.path.basename(fn)
    import tempfile
    dn = tempfile.gettempdir()
    if os.path.isdir(dn):
        fn = os.path.join(dn, fn)
    print(fn)
#    print(outFn.basename, outFn.timestamp)
    #CGSFile.appendHeaderLine(fn, CGSFile.signatures()[-1])

def plotFiles(files, doPlot = True, suffix = None):
    """Expecting a list of tuples."""
    plot = None
    if doPlot and len(files) > 1:
        plot = PlotUncertainty(len(files)+1, suffix = suffix)
    curves = []
    tau = None
    # plot files, sorted by name
    for fn in sorted(files, key = lambda x: x[1]):
        uc = processFitFile(*fn)
        tau = uc[:,0]
        # FIXME: detect/filter curves with same shape but different TAU!
        if len(curves) and curves[-1].shape != tau.shape:
            logging.warn("Ignoring '{}'".format(os.path.join(fn[:1])))
            logging.warn("Dimension mismatch: expected {}, got {}!"
                         .format(curves[-1].shape, tau.shape))
            continue
        curves.append(uc[:,1])
        if plot is not None:
            plot.plot(uc, os.path.basename(fn[0]))
    combined = np.median(np.vstack(curves), axis = 0)
    combined = np.vstack((tau, combined)).T
    minmax = combined[:,1].min(), combined[:,1].max()
    valueRange = minmax[1] - minmax[0]
    if doPlot and plot is not None:
        plot.plot(combined, "median (min: {0:.2e}, range: {1:.2e})"
                            .format(minmax[0], valueRange))
        plot.show()
    # return the normalized summary plot
#    combined[:,1] -= minmax[0]
#    combined[:,1] /= valueRange
    return combined, plot

def fmtAngle(angle):
    return "{0:.0f}".format(angle) + Deg.displayMagnitudeName

def getUncertainties(paths):
    if not isList(paths):
        paths = [paths]
    fitFiles = []
    for dataPath in paths:
        fitFiles += findFitOutput(str(dataPath))
    # sort files by angle first, use orderedDict
    fitFiles.sort(key = lambda x: x[-1])
    filesByAngle = OrderedDict()
    for dn, fn, angle in fitFiles:
        if angle not in filesByAngle:
            filesByAngle[angle] = []
        filesByAngle[angle].append((dn, fn, angle))
    # init the summary plot first
    plots = [PlotUncertainty(len(filesByAngle),
                             suffix = "of median")]
    ucByAngle = OrderedDict()
    for angle, lst in filesByAngle.items():
        combined, separatePlot = plotFiles(lst, doPlot = True,
                suffix = "@" + fmtAngle(angle))
        if separatePlot is not None:
            plots.append(separatePlot)
        # show the angle only, for more than one input
        label = fmtAngle(angle)
        if len(lst) == 1:
            label = os.path.basename(lst[0][0])
        plots[0].plot(combined, label)
        ucByAngle[angle] = combined
    plots[0].show() # finalize the summary plot
    pyplot.show()
    return ucByAngle
    # use ucByAngle for simulation
#    simulate(combined, doPlot)

if __name__ == "__main__":
    assert len(sys.argv) > 1, (
        "Please provide a directory of DLS measurement and McDLS output files!")
    getUncertainties(sys.argv[1:])

# vim: set ts=4 sw=4 sts=4 tw=0:
