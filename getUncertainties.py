#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# retrieves measurement uncertainties from McDLS result files and directories
# which are possibly filtered, cropped and averaged by scattering angle

import sys
import os
import numpy as np
import logging
logging.basicConfig(level = logging.INFO)

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
            goodFiles.append((dirpath, fn))
        fitFiles += goodFiles
    return fitFiles

def processFitFile(dirpath, filename):
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

    @property
    def titleFigure(self): return "Uncertainties"

    @property
    def titleAxes(self): return "extracted " + self.titleFigure.lower()

    @property
    def scalex(self): return "log"

    @property
    def labelx(self): return "tau"

    @property
    def labely(self): return "Uncertainty"

    def __init__(self):
        self._figure = pyplot.figure(figsize = (15, 8), dpi = 80,
                              facecolor = 'w', edgecolor = 'k')
        self._figure.canvas.set_window_title(self.titleFigure)
        self._axes = pyplot.subplot(xscale = self.scalex)
        self._axes.set_title(self.titleAxes)
        self._axes.set_xlabel(self.labelx)
        self._axes.set_ylabel(self.labely)

    def plot(self, uc, label):
        if not isinstance(uc, np.ndarray):
            return
        logging.info("plotting {} {}".format(type(uc), uc.shape))
        self._axes.plot(uc[:,0], uc[:,1], label = label)

    def show(self):
        self._axes.legend(loc = 1, fancybox = True)
        pyplot.setp(pyplot.gca().get_legend().get_texts(), fontsize = 'small')
        PlotResults.plotGrid(self._axes)
        pyplot.show()

def launchPlot(uc):
    from multiprocessing import Process, Queue
    proc = Process(target = plot, args = [uc])
    proc.start()

def simulate(uc):
    from models.dlssphere import DLSSphere
    from dataobj import DLSData
    from utils.units import K, Vis, NM, Deg, Sec
    # set up the dummy DataObj
    dlsData = DLSData(title = "simulated title")
    dlsData.setSampleName("simulated sample")
    dlsData.setFilename("simulated file")
    # usual values from measurements
    dlsData.setTemperature(K.toSi(296))
    dlsData.setViscosity(Vis.toSi(0.932))
    dlsData.setRefractiveIndex(1.33200)
    dlsData.setWavelength(NM.toSi(632.8))
    logging.info("setting up DLSData with:")
    logging.info("    temp: {}, vis: {}, ref.ind.: {}, lambda: {}"
                    .format(dlsData.temperature, dlsData.viscosity,
                            dlsData.refractiveIndex, dlsData.wavelength))
    # init data for a single angle
    dlsData.setAngles(Deg, Deg.toSi(np.array((90.,))))
    logging.info("    angle: {}".format(dlsData.angles))
    dlsData.setTau(Sec, uc[:, 0])
    dlsData.setCorrelation(np.zeros_like(uc[:,1]).reshape((-1, 1)))
    dlsData.initConfig()
    # set up the DLS model
    model = DLSSphere()
    model.radius.setActive(False)
    model.radius.setDisplayValue(50.)
    # calculate finally
    modelData = model.calc(dlsData, np.array((0.,)).reshape((-1, 1)),
                           compensationExponent = 1.)
    dlsData.setCorrelation(modelData.chisqrInt.reshape((-1, 1)), uc[:,1].reshape((-1, 1)))
    from multiprocessing import Process, Queue
    # multithreaded plotting also logs to file
#    pkwargs["logToFile"] = True
#    q = Queue() # allow communication between processes
#    pkwargs["queue"] = q
    # set up the result data, same as in the McSAS class
    result = dict(fitX0 = dlsData.x0.binnedData,
                  fitMeasValMean = np.zeros_like(uc[:,1]).reshape((1, -1)),
                  times = np.array(()))
    proc = Process(target = PlotResults, args = ((result,), dlsData))
    proc.start()

if __name__ == "__main__":
    assert len(sys.argv) > 1, (
        "Please provide a directory of DLS measurement and McDLS output files!")
    dataPath = sys.argv[1]
    fitFiles = findFitOutput(dataPath)
    for fn in fitFiles:
        print(fn)
    plot = PlotUncertainty()
    combined = None
    for fn in fitFiles[0:3]:
        uc = processFitFile(*fn)
        if combined is None:
            combined = uc.copy()
        else:
            combined[:,1] = np.maximum(combined[:,1], uc[:,1])
        plot.plot(uc, fn[-1])
    plot.plot(combined, "combined (maximum)")
    simulate(combined)
    #plot.show()

# vim: set ts=4 sw=4 sts=4 tw=0:
