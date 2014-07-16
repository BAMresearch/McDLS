# -*- coding: utf-8 -*-
# mcsas/plotting.py
"""
Add docstring
"""

import logging
import numpy as np # For arrays
from numpy import size, log10
from cutesnake.utils import isList, isString
import matplotlib
import matplotlib.font_manager as fm
from matplotlib import gridspec

# set up matplotlib.pyplot, do this *before* importing pyplot
try:
    import PySide # verify/test that we have pyside
    # use() gives an error if calling twice
    matplotlib.rcParams['backend'] = 'QT4Agg'
    matplotlib.rcParams['backend.qt4'] = 'PySide'
except ImportError:
    pass # no pyside

from matplotlib.pyplot import (figure, xticks, yticks, errorbar, bar,
        text, plot, grid, legend, title, xlim, ylim, gca, axis,
        close, colorbar, imshow, subplot, axes)
from pylab import show

class PlotResults(object):
    """
    This function plots the output of the Monte-Carlo procedure in two
    windows, with the left window the measured signal versus the fitted
    intensity (on double-log scale), and the righthand window the size
    distribution.
    """

    def __init__(self, allRes, dataset, 
                 axisMargin = 0.3, parameterIdx = None, figureTitle = None,
                 mcsasInstance = None):
        if not isList(allRes) or not len(allRes):
            logging.info("There are no results to plot, breaking up.")
            return

        # set parameters
        self._allRes = allRes
        self._result = allRes[0]
        self._dataset = dataset
        self._axisMargin = axisMargin
        self._parameterIdx = parameterIdx
        self._figureTitle = figureTitle
        self._mcsasInstance = mcsasInstance
        try:
            self._BG = (np.mean(self._result['scalingFactors'][1, :]), 
                np.std(self._result['scalingFactors'][1, :], ddof = 1 ) )
        except:
            self._BG = (0., 0.)

        # set plot font
        fontFamilyArial = ["Arial", "Bitstream Vera Sans", "sans-serif"]
        fontFamilyTimes = ["Times", "DejaVu Serif", "serif"]
        self._plotfont = fm.FontProperties(
                    family = fontFamilyArial)
        self._textfont = fm.FontProperties(
                    family = fontFamilyTimes)

        # set general axes settings:
        self._AxDict = {'axis_bgcolor' : (.95, .95, .95), 
                'xscale' : 'log', 
                'yscale' : 'log',
                }

        # load original Dataset
        self._data = self._dataset.origin
        self._q = self._data[:, 0]
        self._intensity = self._data[:, 1]
        self._intError = self._data[:, 2]

        # number of histograms:
#        self._nHists = mcsasInstance.model.activeParamCount()
        self._nHists = sum((len(p.histograms())
                    for p in mcsasInstance.model.activeParams()))
        self._nR = 1
        # number of ranges: 
        if False and self._nHists > 0: # disabled for testing
            self._ranges = ( mcsasInstance.model.activeParams()[0]
                                .histogram().ranges )
            self._nR = len( self._ranges )
        else:
            self._nR = 1 # no active parameters

        # initialise figure:
        self._fig, self._ah = self.figInit(self._nHists, 
                self._figureTitle, self._nR)

        # show all ranges:
        for rangei in range(self._nR):
            
            #plot intensity fit:
            if dataset.is2d:
                psi = data[:, 3]
                intensity2d = allRes['intensity2d']
                qAxis = ax[self._nHists + nR]
                self.plot2D(self._q, psi, self._intensity, intensity2d, qAxis)

            else:
                # 1D data
                qAxis = self._ah[rangei * 2 * (self._nHists + 1) 
                    + self._nHists + 1]
                fitQ = np.sort(self._result['fitQ'])
                fitIntensity = self._result['fitIntensityMean'][0, 
                        np.argsort(self._result['fitQ'])]
                self.plot1D(self._q, self._intensity, self._intError, 
                        fitQ, fitIntensity, qAxis)

            ## Information on the settings can be shown here:
            InfoAxis = self._ah[rangei * 2 * (self._nHists + 1)]
            # make active:
            self.plotInfo(InfoAxis)
            axes(InfoAxis)
            # axis('tight')

            # plot histograms
            # https://stackoverflow.com/a/952952
            for hi, parHist in enumerate((parHist for histograms in
                        (p.histograms()
                            for p in self._mcsasInstance.model.activeParams())
                        for parHist in histograms)):
                plotPar = parHist.param
                # histogram data:
#                    parStat = parHist.moments
                # get data:
                # histogram axis index:
#                res = self._allRes[parami]
                res = None
                # prep axes:
                hAxis = self._ah[hi + (self._nHists + 1) + 1]
                # plot partial contribution in qAxis
# not yet available:
#                   fitIntensity, fitSTD = parStat.intensity
#                   self.plotPartial(fitQ, fitIntensity, fitSTD, qAxis)
                self.plotHist(res, plotPar, parHist, 
                        hAxis, self._axisMargin, rangei)

                # put the rangeInfo in the plot above
                InfoAxis = self._ah[hi + 1]
                self.plotStats(parHist, self._mcsasInstance, 
                        rangei, self._fig, InfoAxis)

            #plot labels in qAxis:
            axes(qAxis)
            legend(loc = 1, fancybox = True, prop = self._textfont)
        # trigger plot window popup
        self._fig.canvas.draw()
        self._fig.show()

    def plotGrid(self, figh):
        #make axis active:
        axes(figh)
        #plot grid
        grid(lw = 2, color = 'black', alpha = .5, dashes = [1, 6],
             dash_capstyle = 'round', zorder = -1)
            
    def formatRangeInfo(self, parHist, RI, mcsasInstance, weighti = 0):
        """Preformats the rangeInfo results ready for printing"""
        oString = 'Range {} to {}, {}-weighted'.format(
                parHist.lower, parHist.upper,
                parHist.yweight)
        pStat = parHist.moments
        pStatFields = pStat.fields
        pStatFieldNames = pStat.fieldNames()
        for si in np.arange(0,10,2):
            pStatFieldName = pStatFieldNames[si]
            pStatField = pStatFields[si]
            pStatFieldSTD = pStatFields[si + 1]
            oString += '\n {0}:  {1:0.03e} $\pm$ {2:0.03e}'.format(
                    pStatFieldName,
                    pStatField,
                    pStatFieldSTD)

        return oString

    def formatAlgoInfo(self):
        """Preformats the algorithm information ready for printing
        the colons are surrounded by string-marks, to force laTeX rendering"""
        oString = ' Fitting of data$:$ {} '.format(self._figureTitle)
        oString += '\n Q-range$:$ {0:03e} to {1:03e} '.format(
            self._q.min(), self._q.max())
        #oString.append('\n number of datapoints: {}'.format(len(self._q)))
        oString += '\n Active parameters$:$ {}, ranges: {} '.format(
            self._mcsasInstance.model.activeParamCount(), self._nR)
        oString += '\n Background level: {0:03e} $\pm$ {1:03e}'.format(
                self._BG[0], self._BG[1])

        return oString

    def setAxis(self, ah):
        # self.setAxis font and ticks
        ah.set_yticklabels(ah.get_yticks(), fontproperties = self._plotfont,
                size = 'large')
        ah.set_xticklabels(ah.get_xticks(), fontproperties = self._plotfont,
                size = 'large')
        ah.set_xlabel(ah.get_xlabel(), fontproperties = self._textfont,
                size = 'x-large')
        ah.set_ylabel(ah.get_ylabel(), fontproperties = self._textfont,
                size = 'x-large')
        ah.spines['bottom'].set_lw(2)
        ah.spines['top'].set_lw(2)
        ah.spines['left'].set_lw(2)
        ah.spines['right'].set_lw(2)
        ah.tick_params(axis = 'both', colors = 'black', width = 2,
                which = 'major', direction = 'in', length = 6)
        ah.tick_params(axis = 'x', colors = 'black', width = 2,
                which = 'minor', direction = 'in', length = 3)
        ah.tick_params(axis = 'y', colors = 'black', width = 2,
                which = 'minor', direction = 'in', length = 3)
        locs, labels = xticks()
        xticks(locs, map(lambda x: "%g" % x, locs))
        locs, labels = yticks()
        yticks(locs, map(lambda x: "%g" % x, locs))
        return ah

    def figInit(self, nHists, figureTitle, nR = 1):
        """initialize figure and initialise axes using GridSpec.
        Each rangeinfo (nR) contains two rows and nHists + 1 columns.
        the top row axes are for placing text objects: settings and stats.
        The bottom row axes are for plotting the fits and the histograms
        TODO: add settings to window title? (next to figure_xy)"""
        ah = list() #list of axes handles from top left to bottom right.

        fig = figure(figsize = (7 * (nHists + 1), 7 * nR), dpi = 80,
                     facecolor = 'w', edgecolor = 'k')
        if isString(figureTitle):
            fig.canvas.set_window_title(figureTitle)

        gs = gridspec.GridSpec(2 * nR, nHists + 1,
                height_ratios = np.tile([1,6],nR ) )
        #update margins
        gs.update(left = 0.08, bottom = 0.10,
                            right = 0.96, top = 0.95,
                            wspace = 0.23, hspace = 0.13)

        for ai in range((nHists + 1) * nR * 2 ):
            #initialise axes 
            ah.append(subplot(gs[ai]))
            if ai%((nHists + 1) * 2) < (nHists + 1) : 
                #text box settings:
                textAxDict = {
                        'frame_on' : False,
                        'yticks' : [],
                        'xticks' : [],
                        'ylim' : [0., 1.],
                        'xlim' : [0., 1.],
                        }
                ah[-1].update(textAxDict)
        return fig, ah

    def plot2D(self, q, psi, intensity, intensity2d, qAxis):
        """plots 2D data and fit"""
        # 2D data
        # we need to recalculate the result in two dimensions
        intShow = intensity.copy()
        # quadrant 1 and 4 are simulated data, 2 and 3 are measured data
        intShow[(psi >   0) * (psi <=  90)] = intensity2d[
                (psi >   0) * (psi <=  90)]
        intShow[(psi > 180) * (psi <= 270)] = intensity2d[
                (psi > 180) * (psi <= 270)]
        xmidi = int(round(size(q, 1)/2))
        ymidi = int(round(size(q, 0)/2))
        QX = np.array([-q[ymidi, 0], q[ymidi, -1]])
        QY = np.array([-q[0, xmidi], q[-1, xmidi]])
        extent = (QX[0], QX[1], QY[0], QY[1])

        # indexing probably wrong:
        qAxis.update( axisbg = (.95, .95, .95),
                               xlim = QX, ylim = QY, xlabel = 'q_x, $\times 10^9$ 1/m',
                               ylabel = 'q_y, $\times 10^9$  1/m')
        imshow(log10(intShow), extent = extent, origin = 'lower')
        qAxis = self.setAxis(qAxis)
        colorbar()
        title('Measured vs. Fitted intensity',
              fontproperties = self._textfont, size = 'large')
        # reapply limits, necessary for some reason:
        xlim(QX)
        ylim(QY)

    def plotPartial(self, fitQ, fitIntensity, fitSTD, qAxis):
        """plots 1D data and fit"""
        #make active:
        axes(qAxis)
        plot(fitQ, fitIntensity, 'b-', lw = 1, label = 'MC partial intensity')

    def plot1D(self, q, intensity, intError, fitQ, fitIntensity, qAxis):
        #settings for Q-axes (override previous settings where appropriate):
        xLim = (q.min() * (1 - self._axisMargin), 
                q.max() * (1 + self._axisMargin))
        yLim = (intensity[intensity != 0].min() * (1 - self._axisMargin), 
                intensity.max() * (1 + self._axisMargin))
        qAxDict = self._AxDict.copy()
        qAxDict.update({
                'xlim' : xLim,
                'ylim' : yLim,
                'xlabel' : 'q, $\\times 10^9 m^{-1}$', 
                'ylabel' : 'intensity, $(m\, sr)^{-1}$'
                })

        """plots 1D data and fit"""
        #make active:
        axes(qAxis)
        qAxis.update(qAxDict)
        qAxis = self.setAxis(qAxis)
        errorbar(q, intensity, intError, zorder = 2, fmt = 'k.',
                 ecolor = 'k', elinewidth = 2, capsize = 4, ms = 5,
                 label = 'Measured intensity', lw = 2,
                 solid_capstyle = 'round', solid_joinstyle = 'miter')
        self.plotGrid(qAxis)
        plot(fitQ, fitIntensity, 'r-', lw = 3, 
                label = 'MC Fit intensity', zorder = 4)
        try:
            plot(fitQ, self._BG[0] + 0*fitQ,
                 'g-', linewidth = 3,
                 label = 'MC Background level:\n\t ({0:03.3g})' .format(
                     self._BG[0]), zorder = 3)
        except:
            logging.error('could not plot background')
            pass
        title('Measured vs. Fitted intensity',
              fontproperties = self._textfont, size = 'large')
        # reapply limits, necessary for some reason:
        xlim(xLim)
        ylim(yLim)


    def plotInfo(self, InfoAxis):
        """plots the range statistics in the small info axes above plots"""
        delta = 0.001 #minor offset
        #make active:
        axes(InfoAxis)
        #show volume-weighted info:
        ovString = self.formatAlgoInfo()
        tvObj = text(0. - delta, 0. + delta, ovString,
                family = "monospace", size = "small", 
                horizontalalignment = 'center',
                multialignment = 'center',
                verticalalignment = 'center')
        self._fig.show()
        axis('tight')

    def plotStats(self, parHist, mcsasInstance, rangei, fig, InfoAxis):
        """plots the range statistics in the small info axes above plots"""
        delta = 0.001 #minor offset
        #make active:
        axes(InfoAxis)
        #show volume-weighted info:
        ovString = self.formatRangeInfo(parHist, rangei, mcsasInstance, weighti = 0)
        tvObj = text(0. - delta, 0. + delta, ovString, bbox = 
                {'facecolor' : 'white', 'alpha': 0.95},
                family = "monospace", size = "small", 
                horizontalalignment = 'right',
                multialignment = 'right',
                verticalalignment = 'center')
        fig.show()
        axis('tight')
        return # stop here
        #get bounding box
        #bb = tvObj.get_window_extent()
        #width = bb.width
        #height = bb.height
        #print('width: {}, height: {}'.format(width, height))
        aLim = axis()
        #axis('tight')
        #add number-weighted info:
        onString = self.formatRangeInfo(parHist, rangei, mcsasInstance, weighti = 1)
        tnObj = text(0. + delta, 0. + delta, onString, bbox = 
                {'facecolor' : 'white', 'alpha': 0.95},
                family = "monospace", size = "small", 
                horizontalalignment = 'left',
                multialignment = 'right',
                verticalalignment = 'center')
        fig.show()
        axis('tight')


    def plotHist(self, res, plotPar, parHist, hAxis, axisMargin, rangei):
        """histogram plot"""
        #make active:
        axes(hAxis)
#        hRange = parHist.ranges[rangei]

        histXLowerEdge = parHist.xLowerEdge
        histXMean = parHist.xMean
        histXWidth = parHist.xWidth
        # plot volume weighted by default, both would be good
        # can we plot both weightings? perhaps, with different colors?
        # e.g. red/orange (current) and blue/lightblue?
        volHistYMean = parHist.bins.mean
        volHistMinReq = parHist.observability
        volHistYStd = parHist.bins.std
        #elif params[parameterId[parami]].histogram().hasWeighting('num'):
        #    volHistYMean = res['numberHistogramYMean']
        #    volHistMinReq = res['numberHistogramMinimumRequired']
        #    volHistYStd = res['numberHistogramYStd']
        #else: 
        #    "Incorrect value for histWeighting: "\
        #    "should be either 'volume' or 'number'"

        #get information for labels:
        plotTitle = plotPar.displayName()
        xLabel = '{}, {}'.format(plotPar.name(), plotPar.suffix())

        if parHist.xscale == 'log':
            xLim = (histXLowerEdge.min() * (1 - self._axisMargin), 
                    histXLowerEdge.max() * (1 + self._axisMargin))
            xScale = 'log'
        else:
            xLim = (histXLowerEdge.min() - (1 - self._axisMargin)
                    * histXLowerEdge.min(), 
                    histXLowerEdge.max() * (1 + self._axisMargin))
            xScale = 'linear'

        yLim = (0, volHistYMean.max() * (1 + self._axisMargin) )
        # histogram axes settings:
        hAxDict = self._AxDict.copy()
        # change axis settings:
        hAxDict.update({
            'xlim' : xLim,
            'ylim' : yLim,
            'xlabel' : xLabel,
            'xscale' : xScale,
            'yscale' : 'linear',
            'ylabel' : '[Rel.] Volume Fraction' })
        # update axes settings:
        hAxis.update(hAxDict)
        # change axis settings not addressible through dictionary:
        hAxis = self.setAxis(hAxis)
        #plot grid
        self.plotGrid(hAxis)

        # fill axes
        # plot inactive histogram:
        bar(histXLowerEdge[0:-1], volHistYMean, 
                width = histXWidth, color = 'grey',
                edgecolor = 'black', linewidth = 0.5, zorder = 2, alpha = 0.5,
                )
        # plot active histogram:
        validi = (histXLowerEdge >= parHist.lower) * (histXLowerEdge <= parHist.upper)
        validi[-1] = 0
        if not (validi.sum()==0):
            bar(histXLowerEdge[validi], volHistYMean[validi[0:-1]], 
                    width = histXWidth[validi[0:-1]], color = 'orange',
                    edgecolor = 'black', linewidth = 1, zorder = 2,
                    label = 'MC size histogram')
        # plot observability limit
        plot(histXMean, volHistMinReq, 'ro', 
                ms = 5, markeredgecolor = 'r',
                label = 'Minimum visibility limit', zorder = 3)
        # plot inactive uncertainties
        errorbar(histXMean, volHistYMean, volHistYStd,
                zorder = 4, fmt = 'k.', ecolor = 'k',
                elinewidth = 0.5, capsize = 3, ms = 0, lw = 0.5,
                solid_capstyle = 'round', solid_joinstyle = 'miter')
        # plot active uncertainties
        errorbar(histXMean[validi[0:-1]], volHistYMean[validi[0:-1]], 
            volHistYStd[validi[0:-1]],
                zorder = 4, fmt = 'k.', ecolor = 'k',
                elinewidth = 2, capsize = 4, ms = 0, lw = 2,
                solid_capstyle = 'round', solid_joinstyle = 'miter')
        legend(loc = 1, fancybox = True, prop = self._textfont)
        title(plotTitle, fontproperties = self._textfont,
              size = 'large')
        xlim(xLim)



# vim: set ts=4 sts=4 sw=4 tw=0:
