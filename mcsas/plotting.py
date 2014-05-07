# -*- coding: utf-8 -*-
# mcsas/plotting.py
"""
Add docstring
"""

import numpy as np # For arrays
from numpy import size, log10
from cutesnake.utils import isList, isString
import matplotlib
import matplotlib.font_manager as fm
from matplotlib import gridspec
from matplotlib.pyplot import (figure, xticks, yticks, errorbar, bar, 
        text, plot, grid, legend, title, xlim, ylim, gca, axis,
        close, colorbar, imshow, subplot, axes)
from pylab import show

def plotResults(allRes, dataset, 
                axisMargin = 0.3, parameterIdx = None, figureTitle = None,
                mcsasInstance = None):
    """
    This function plots the output of the Monte-Carlo procedure in two
    windows, with the left window the measured signal versus the fitted
    intensity (on double-log scale), and the righthand window the size
    distribution.
    """

    try:
        import PySide # verify/test that we have pyside
        # use() gives an error if calling twice
        matplotlib.rcParams['backend'] = 'QT4Agg'
        matplotlib.rcParams['backend.qt4'] = 'PySide'
    except ImportError:
        pass # no pyside

    def formatRangeInfo(parameter, RI, mcsasInstance, weighti = 0):
        """Preformats the rangeInfo results ready for printing"""
        weightings = parameter.weighting()
        weighting = weightings[weighti]
        oString = 'Range {} to {}, {}-weighted'.format(
                parameter.ranges[RI][0],
                parameter.ranges[RI][1],
                weighting)
        pStat = parameter.stats[RI][weighti]
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

        print('{}'.format(oString))
        return oString

    def setAxis(ah):
        """Sets the axes Parameters."""
        plotfont = fm.FontProperties(
                    family = fontFamilyArial)
        textfont = fm.FontProperties(
                    family = fontFamilyTimes)
        # setAxis font and ticks
        ah.set_yticklabels(ah.get_yticks(), fontproperties = plotfont,
                size = 'large')
        ah.set_xticklabels(ah.get_xticks(), fontproperties = plotfont,
                size = 'large')
        ah.set_xlabel(ah.get_xlabel(), fontproperties = textfont,
                size = 'x-large')
        ah.set_ylabel(ah.get_ylabel(), fontproperties = textfont,
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

    def figInit(nHists, figureTitle, nR = 1):
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

    def plot2D(q, psi, intensity, intensity2d, qAxis):
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
                               xlim = QX, ylim = QY, xlabel = 'q_x, 1/m',
                               ylabel = 'q_y, 1_m')
        imshow(log10(intShow), extent = extent, origin = 'lower')
        qAxis = setAxis(qAxis)
        colorbar()

    def plot1D(q, intensity, intError, fitQ, fitIntensity, qAxis):
        """plots 1D data and fit"""
        #make active:
        axes(qAxis)
        qAxis.update(qAxDict)
        qAxis = setAxis(qAxis)
        errorbar(q, intensity, intError, zorder = 2, fmt = 'k.',
                 ecolor = 'k', elinewidth = 2, capsize = 4, ms = 5,
                 label = 'Measured intensity', lw = 2,
                 solid_capstyle = 'round', solid_joinstyle = 'miter')
        grid(lw = 2, color = 'black', alpha = .5, dashes = [1, 6],
             dash_capstyle = 'round', zorder = -1)
        plot(fitQ, fitIntensity, 'r-', lw = 3, label = 'MC Fit intensity', zorder = 4)
        try:
            plot(aq, np.mean(result['scalingFactors'][1, :]) + 0*aq,
                 'g-', linewidth = 3,
                 label = 'MC Background level:\n\t ({0:03.3g})'
                         .format(np.mean(result['scalingFactors'][1, :])),
                 zorder = 3)
        except:
            pass
        legend(loc = 1, fancybox = True, prop = textfont)

    def plotStats(parHist, mcsasInstance, rangei, fig, InfoAxis):
        """plots the range statistics in the small info axes above plots"""
        delta = 0.001 #minor offset
        #make active:
        axes(InfoAxis)
        #show volume-weighted info:
        ovString = formatRangeInfo(parHist, rangei, mcsasInstance, weighti = 0)
        tvObj = text(0. - delta, 0. + delta, ovString, bbox = 
                {'facecolor' : 'white', 'alpha': 0.5},
                family = "monospace", size = "small", 
                horizontalalignment = 'right',
                multialignment = 'right',
                verticalalignment = 'center')
        fig.show()
        #get bounding box
        #bb = tvObj.get_window_extent()
        #width = bb.width
        #height = bb.height
        #print('width: {}, height: {}'.format(width, height))
        aLim = axis()
        #axis('tight')
        #add number-weighted info:
        onString = formatRangeInfo(parHist, rangei, mcsasInstance, weighti = 1)
        tnObj = text(0. + delta, 0. + delta, onString, bbox = 
                {'facecolor' : 'white', 'alpha': 0.5},
                family = "monospace", size = "small", 
                horizontalalignment = 'left',
                multialignment = 'right',
                verticalalignment = 'center')
        fig.show()
        axis('tight')


    def plotHist(res, plotPar, parHist, hAxis, axisMargin):
        """histogram plot"""
        #make active:
        axes(hAxis)

        histXLowerEdge = res['histogramXLowerEdge']
        histXMean = res['histogramXMean']
        histXWidth = res['histogramXWidth']
        # plot volume weighted by default, both would be good
        # can we plot both weightings? perhaps, with different colors?
        # e.g. red/orange (current) and blue/lightblue?
        volHistYMean = res['volumeHistogramYMean']
        volHistMinReq = res['volumeHistogramMinimumRequired']
        volHistYStd = res['volumeHistogramYStd']
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

        if parHist.scaleX == 'log':
            xLim = (histXLowerEdge.min() * (1 - axisMargin), 
                    histXLowerEdge.max() * (1 + axisMargin))
            xScale = 'log'
        else:
            xLim = (histXLowerEdge.min() - (1 - axisMargin)
                    * histXLowerEdge.min(), 
                    histXLowerEdge.max() * (1 + axisMargin))
            xScale = 'linear'

        yLim = (0, volHistYMean.max() * (1 + axisMargin) )
        # change axis settigns:
        hAxDict.update({
            'xlim' : xLim,
            'ylim' : yLim,
            'xlabel' : xLabel,
            'xscale' : xScale,
            'yscale' : 'linear',
            'ylabel' : '[Rel.] Volume Fraction' })
        # update axes settings:
        hAxis.update(hAxDict)
        # store in list of histogram axes (maybe depreciated soon):
        sizeAxis.append(hAxis)
        # change axis settings not addressible through dictionary:
        sizeAxis[parami] = setAxis(sizeAxis[parami])
        # fill axes
        # plot histogram:
        bar(histXLowerEdge[0:-1], volHistYMean, 
                width = histXWidth, color = 'orange',
                edgecolor = 'black', linewidth = 1, zorder = 2,
                label = 'MC size histogram')
        # plot observability limit
        plot(histXMean, volHistMinReq, 'ro', 
                ms = 5, markeredgecolor = 'r',
                label = 'Minimum visibility limit', zorder = 3)
        # plot uncertainties
        errorbar(histXMean, volHistYMean, volHistYStd,
                zorder = 4, fmt = 'k.', ecolor = 'k',
                elinewidth = 2, capsize = 4, ms = 0, lw = 2,
                solid_capstyle = 'round', solid_joinstyle = 'miter')
        legend(loc = 1, fancybox = True, prop = textfont)
        title(plotTitle, fontproperties = textfont,
              size = 'large')
        # reapply limits in x
        xlim((histXLowerEdge.min() * (1 - axisMargin),
              histXLowerEdge.max() * (1 + axisMargin)))


    # set plot font
    fontFamilyArial = ["Arial", "Bitstream Vera Sans", "sans-serif"]
    fontFamilyTimes = ["Times", "DejaVu Serif", "serif"]
    plotfont = fm.FontProperties(
                #size = 'large',
                family = fontFamilyArial)
    textfont = fm.FontProperties(
                #size = 'large',
                family = fontFamilyTimes)

    # load original Dataset
    data = dataset.origin
    q = data[:, 0]
    intensity = data[:, 1]
    intError = data[:, 2]

    if not isList(allRes) or not len(allRes):
        logging.info("There are no results to plot, breaking up.")
        return
    result = allRes[0]

    # check how many result plots we need to generate, and find the 
    # indices to the to-plot paramters
    # number of histograms:
    nHists = mcsasInstance.model.activeParamCount()
    # number of ranges: 
    if nHists > 0:
        nR = len( mcsasInstance.model.activeParams()[0].histogram().ranges )
    else:
        nR = 1 # no active parameters
    # initialise figure:
    fig, ah = figInit(nHists, figureTitle, nR)

    #general axes settings:
    AxDict = {'axis_bgcolor' : (.95, .95, .95), 
            'xscale' : 'log', 
            'yscale' : 'log',
            }

    #settings for Q-axes (override previous settings where appropriate):
    xLim = (q.min() * (1 - axisMargin), q.max() * (1 + axisMargin))
    yLim = (intensity[intensity != 0].min() * (1 - axisMargin), 
            intensity.max() * (1 + axisMargin))
    qAxDict = AxDict.copy()
    qAxDict.update({
            'xlim' : xLim,
            'ylim' : yLim,
            'xlabel' : 'q, 1/m', 
            'ylabel' : 'intensity, 1/(m sr)'
            })

    # quickFix for now, to be modified later for showing more ranges:
    rangei = 0
    
    #plot intensity fit:
    if dataset.is2d:
        psi = data[:, 3]
        intensity2d = allRes['intensity2d']
        qAxis = ax[nHists + nR]
        plot2D(q, psi, intensity, intensity2d, qAxis)

    else:
        # 1D data
        qAxis = ah[(rangei + 1) * (nHists + 1) ]
        fitQ = np.sort(result['fitQ'])
        fitIntensity = result['fitIntensityMean'][0, 
                np.argsort(result['fitQ'])]
        plot1D(q, intensity, intError, fitQ, fitIntensity, qAxis)

    title('Measured vs. Fitted intensity',
          fontproperties = textfont, size = 'large')
    # reapply limits, necessary for some reason:
    xlim(xLim)
    ylim(yLim)

    # Information on the settings can be shown here:
    InfoAxis = ah[rangei * (nHists + 1)]
    # make active:
    axes(InfoAxis)
    text(0,0,'Algorithm settings \n go here...', bbox = 
            {'facecolor' : 'white', 'alpha': 0.5},
            fontproperties = textfont)
    axis('tight')

    sizeAxis = list()
    # plot histograms
    # histogram axes settings:
    hAxDict = AxDict.copy()
    for parami, plotPar in enumerate(mcsasInstance.model.activeParams()):
        # histogram data:
        parHist = plotPar.histogram()
        # get data:
        # histogram axis index:
        res = allRes[parami]
        # prep axes:
        hAxis = ah[(rangei + 1) * nHists + 2 + parami]

        plotHist(res, plotPar, parHist, hAxis, axisMargin)

        #put the rangeInfo in the plot above
        InfoAxis = ah[(rangei) * nHists + 1 + parami]
        plotStats(parHist, mcsasInstance, rangei, fig, InfoAxis)

    # trigger plot window popup
    fig.show()

# vim: set ts=4 sts=4 sw=4 tw=0:
