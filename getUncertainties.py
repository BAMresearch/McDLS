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
import tempfile
from collections import OrderedDict
import numpy as np
import logging
logging.basicConfig(level = logging.INFO)
#np.seterr(all='raise', under = 'ignore')
from multiprocessing import Process, Pipe
import matplotlib
import matplotlib.pyplot as pyplot

from utils import isString, isList, mcopen
from utils.units import Deg, NM3
from mcsas.plotting import PlotResults
from log import timestampFormatted

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
    with mcopen(filepath, 'rb') as fd:
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

class PlotUncertainty(object):
    _figure = None
    _axes = None
    _cmap = None # colormap
    _idx = None  # number of plots, counting
    _prefix = None # title prefix
    _suffix = None # title suffix

    @property
    def titleFigure(self):
        title = "Uncertainties"
        if isString(self._suffix):
            title = "{} {}".format(title, self._suffix)
        if isString(self._prefix):
            title = "{} {}".format(self._prefix, title)
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
        if isString(self._prefix):
            lbl = "{} {}".format(self._prefix, lbl)
        return lbl

    def __init__(self, count = 0, prefix = None, suffix = None, conn = None):
        self._prefix = prefix
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
        while True:
            # waiting for commands to call plot() or show()
            cmd, args, kwargs = conn.recv()
            func = getattr(self, cmd, None)
            if func is not None:
                func(*args, **kwargs)
                # show() stops the waiting here,
                # pyplot in show() blocks itself (win&lin)
                if "show" in cmd: break # done here

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
        pyplot.show()

class PlotUncertaintyProcess(Process):
    _parentConn = None
    _childConn  = None

    def __init__(self, *args, **kwargs):
        self._parentConn, self._childConn = Pipe()
        kwargs["conn"] = self._childConn
        super(PlotUncertaintyProcess, self).__init__(target = PlotUncertainty,
                                                     kwargs = kwargs)
        self.start()

    def plot(self, *args, **kwargs):
        self._parentConn.send(("plot", args, kwargs))

    def show(self, *args, **kwargs):
        self._parentConn.send(("show", args, kwargs))

def volume(pop):
    """Returns the volume of a given list of radii (aka population)"""
    def vol(r):
        return 3./4. * np.pi * r**3
    return sum(vol(r) for r in pop)

def simulate(uc, label, doPlot = True):
    try:
        assert len(uc.shape) == 2
    except (AttributeError, AssertionError):
        return None
    from models.dlssphere import DLSSphere
    from dataobj import DLSData
    from utils.units import K, Vis, NM, Deg, Sec, MSec
    # set up the dummy DataObj
    dlsData = DLSData(title = "simulated")
    dlsData.setSampleDescription("{} uncertainty".format(label))
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
#    lst = [NM.toSi(x) for x in (13., 15., 15., 17.)] # RM8012, 15nm radius
    volRatio = 1./20 # desired volume ratio
    populations = []
    # FD102, 20nm
    populations.append([NM.toSi(x) for x in ( 9., 10., 10., 11.)])
    # FD102, 80nm
    populations.append([NM.toSi(x) for x in (36., 40., 40., 44.)])
    for i, pop in enumerate(populations):
        logging.info("simulated population {}: {} ({})"
                        .format(i, NM.toDisplay(pop), NM.displayMagnitudeName))
    # calculate amount of each population to reach the desired ratio
    volumes = [volume(pop) for pop in populations]
    # assumes the small population first, the larger ones next
    logging.info("population volumes: {} ({})".format(NM3.toDisplay(volumes),
                                                      NM3.displayMagnitudeName))
    if len(volumes) > 1:
        factorSmall = (1./volRatio -1.) * volumes[1]/volumes[0]
        logging.info("Adjusting smaller population count by factor {:.4g} "
                     "to reach the desired volume ratio of {:.4g}:{}."
                     .format(factorSmall, (1./volRatio -1.), "1"))
        # increase numbers of small population
        decimals = factorSmall%1
        tol = 1e-2
        if decimals < tol or (1.-decimals) < tol:
            populations[0] *= int(factorSmall)
        else: # multiply all count by inverse tolerance
            populations[0] *= int(factorSmall/tol)
            populations[1] *= int(1./tol)
        logging.info("Overall populations count: {}"
                     .format(sum([len(pop) for pop in populations])))
    volumes = []
    for i, pop in enumerate(populations):
        volumes.append(volume(pop))
        logging.info("new volume of population {}: {} {}".format(
            i, NM3.toDisplay(volumes[-1]), NM3.displayMagnitudeName))
    logging.info("volume ratio: {}".format(
        ":".join(["{0:.4g}".format(v) for v in (np.array(volumes) / min(volumes))])))
    contrib = np.array([r for pop in populations for r in pop]).reshape((-1, 1))
    storePopulation(contrib)
    modelData = model.calc(dlsData, contrib, compensationExponent = 1.)
    corr = modelData.chisqrInt.reshape((-1, 1))
#    corrU = uc[:,1].reshape((-1, 1))               # uncertainty from data
#    corrU = np.ones_like(uc[:,1]).reshape((-1, 1)) # no uncertainty, =1
    corrU = None
    dlsData.setCorrelation(corr, corrU)
    return dlsData

def storePopulation(pop):
    fn = os.path.join(tempfile.gettempdir(),
                      "sim_population_{}".format(timestampFormatted()))
    level = 0.2
    if level > 0.:
        err = np.random.uniform(-level, level, pop.shape) + 1.
        logging.info("Storing population with an error of {:g}% "
                     "by multiplying with values in [{:g},{:g}]."
                     .format(level*100., err.min(), err.max()))
        np.savetxt(fn, pop * err)
    else:
        np.savetxt(fn, pop)
    logging.info("Simulated population stored to 'file://{}'.".format(fn))

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

from datafile import AsciiFile
from gui.calc import OutputFilename

class DummyDataSet(object):
    """Just for proper output file formatting."""
    title = u"uncertainties"

def storePlotData(fileKey, dataDict):
    outFn = OutputFilename(DummyDataSet, createDir = False)
    fn = outFn.filenameVerbose(fileKey, fileKey, extension = ".dat")
    columnNames = list(dataDict.viewkeys())
    # write header:
    AsciiFile.writeHeaderLine(fn, [' "{}"'.format(cn) for cn in columnNames])
    # (additional header lines can be appended if necessary)
    data = np.vstack([dataDict[cn] for cn in columnNames]).T
    # append to the header, do not overwrite:
    AsciiFile.appendFile(fn, data)

def plotFiles(files, doPlot = True, suffix = None):
    """Expecting a list of tuples."""
    plot = None
    if doPlot and len(files) > 1:
        plot = PlotUncertaintyProcess(len(files)+1, suffix = suffix)
    dataDict = OrderedDict()
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
        lbl = os.path.basename(fn[0])
        if not len(dataDict):
            dataDict["tau"] = tau
        dataDict[lbl] = curves[-1]
        if plot is not None:
            plot.plot(uc, lbl)
    combined = np.median(np.vstack(curves), axis = 0)
    combined = np.vstack((tau, combined)).T
    minmax = combined[:,1].min(), combined[:,1].max()
    valueRange = minmax[1] - minmax[0]
    lbl = ("median (min: {0:.2e}, max: {0:.2e}, range: {1:.2e})"
           .format(minmax[0], minmax[1], valueRange))
    dataDict[lbl] = combined[:,1] # append the median
    if doPlot and plot is not None:
        plot.plot(combined, lbl)
        plot.show()
    # store this data
    storePlotData(suffix, dataDict)
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
    if not len(fitFiles):
        logging.error("No results found in path '{}'!".format(paths))
        return None
    # sort files by angle first, use orderedDict
    fitFiles.sort(key = lambda x: x[-1])
    filesByAngle = OrderedDict()
    for dn, fn, angle in fitFiles:
        if angle not in filesByAngle:
            filesByAngle[angle] = []
        filesByAngle[angle].append((dn, fn, angle))
    # init the summary plot first
    titlePrefix = "median of"
    plots = [PlotUncertaintyProcess(len(filesByAngle),
                                    prefix = titlePrefix)]
    dataDict = OrderedDict()  # for file output
    ucByAngle = OrderedDict() # for use in data simulation
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
        ucByAngle[angle] = combined, label
        if not len(dataDict):
            dataDict["tau"] = combined[:, 0]
        dataDict[label] = combined[:, 1]
    if "median" in titlePrefix.lower():
        titlePrefix = titlePrefix.split()[0]
    storePlotData(titlePrefix, dataDict)
    plots[0].show() # finalize the summary plot
    return ucByAngle

if __name__ == "__main__":
    assert len(sys.argv) > 1, (
        "Please provide a directory of DLS measurement and McDLS output files!")
    getUncertainties(sys.argv[1:])

# vim: set ts=4 sw=4 sts=4 tw=0:
