#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# retrieves measurement uncertainties from McDLS result files and directories
# which are possibly filtered, cropped and averaged by scattering angle

import sys
import os

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
    for fn in fitFiles:
        print(fn)
    return fitFiles

import numpy as np

def processFitFile(dirpath, filename):
    print("processFitFile", dirpath, filename)
    filepath = os.path.join(dirpath, filename)
    lines = None
    with open(filepath, 'rb') as fd:
        lines = [l.decode("utf8") for l in fd.readlines()]
    if lines is None or len(lines) < 2:
        return
    fields = lines[0].split()
    colidx = [i for i, field in enumerate(fields) if "unbinnedStd" in field]
    if not len(colidx) or not "fitX0" in fields[0]:
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
from matplotlib.pyplot import (figure, subplot)

def plot(uc):
    if not isinstance(uc, np.ndarray):
        return
    print("plot", type(uc), uc.shape)
    fig = figure(figsize = (7, 7), dpi = 80,
                 facecolor = 'w', edgecolor = 'k')
    axes = subplot(xscale = 'log')
    axes.plot(uc[:,0], uc[:,1])
#    fig.canvas.draw()
#    fig.show()
    matplotlib.pyplot.show()

if __name__ == "__main__":
    assert len(sys.argv) > 1, (
        "Please provide a directory of DLS measurement and McDLS output files!")
    dataPath = sys.argv[1]
    fitFiles = findFitOutput(dataPath)
    for fn in fitFiles[0:3]:
        uc = processFitFile(*fn)
        plot(uc)

# vim: set ts=4 sw=4 sts=4 tw=0:
