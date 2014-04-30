# -*- coding: utf-8 -*-
# mcsas/plotting.py
"""
Add docstring
"""

import numpy # For arrays
from numpy import size, log10
from cutesnake.utils import isList, isString
from utils.parameter import Parameter
#from mcsasparameters import McSASParameters 

def plotResults(allRes, dataset, params,
                axisMargin = 0.3, parameterIdx = None, figureTitle = None,
                mcsasInstance = None):
    """
    This function plots the output of the Monte-Carlo procedure in two
    windows, with the left window the measured signal versus the fitted
    intensity (on double-log scale), and the righthand window the size
    distribution.
    """

    import matplotlib
    try:
        import PySide # verify/test that we have pyside
        # use() gives an error if calling twice
        matplotlib.rcParams['backend'] = 'QT4Agg'
        matplotlib.rcParams['backend.qt4'] = 'PySide'
    except ImportError:
        pass # no pyside
    import matplotlib.font_manager as fm
    from matplotlib.pyplot import (figure, xticks, yticks, errorbar, bar,
                                   plot, grid, legend, title, xlim, gca,
                                   close, colorbar, imshow)
    from pylab import show

    def setAxis(ah):
        """Sets the axes Parameters."""
        plotfont = fm.FontProperties(
                    # this only works for macs, doesn't it?
                    # family = 'Courier New Bold',
                    # fname = '/Library/Fonts/Courier New Bold.ttf')
                    family = fontFamilyArial)
        textfont = fm.FontProperties(
                    # Baskerville.ttc does not work when saving to eps
                    # family = 'Times New Roman',
                    # fname = '/Library/Fonts/Times New Roman.ttf')
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

    def figInit(nHists, figureTitle):
        # initialize figure 
        # TODO: add settings to window title? (next to figure_xy)
        fig = figure(figsize = (7*(nHists+1), 7), dpi = 80,
                     facecolor = 'w', edgecolor = 'k')
        if isString(figureTitle):
            fig.canvas.set_window_title(figureTitle)
        return fig

    # set plot font
    fontFamilyArial = ["Arial", "Bitstream Vera Sans", "sans-serif"]
    fontFamilyTimes = ["Times", "DejaVu Serif", "serif"]
    plotfont = fm.FontProperties(
                size = 'large',
                family = fontFamilyArial)
    textfont = fm.FontProperties(
                size = 'large',
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
    nHists = mcsasInstance.model.activeParamCount()

    params = [Parameter(**attr) for attr in params]
    if parameterIdx is None: # use 'is': None is a singleton in python
        parameterId = [i for i, p in enumerate(params) if p.isActive()]
    else:
        parameterId = [parameterIdx]

    fig = figInit(nHists, figureTitle):

    #plot intensity fit:
    if dataset.is2d:
        # 2D data
        psi = data[:, 3]
        # we need to recalculate the result in two dimensions
        # done by gen2DIntensity function
        intensity2d = allRes['intensity2d']
        intShow = intensity.copy()
        # quadrant 1 and 4 are simulated data, 2 and 3 are measured data
        intShow[(psi >   0) * (psi <=  90)] = intensity2d[
                (psi >   0) * (psi <=  90)]
        intShow[(psi > 180) * (psi <= 270)] = intensity2d[
                (psi > 180) * (psi <= 270)]
        xmidi = int(round(size(q, 1)/2))
        ymidi = int(round(size(q, 0)/2))
        QX = numpy.array([-q[ymidi, 0], q[ymidi, -1]])
        QY = numpy.array([-q[0, xmidi], q[-1, xmidi]])
        extent = (QX[0], QX[1], QY[0], QY[1])

        qAxis = fig.add_subplot(1, (nHists+1), 1, axisbg = (.95, .95, .95),
                               xlim = QX, ylim = QY, xlabel = 'q_x, 1/m',
                               ylabel = 'q_y, 1_m')
        imshow(log10(intShow), extent = extent, origin = 'lower')
        qAxis = setAxis(qAxis)
        colorbar()
    else:
        #1D data
        qAxis = fig.add_subplot(
                    1, (nHists+1), 1,
                    axisbg = (.95, .95, .95),
                    xlim = (q.min() * (1 - axisMargin),
                            q.max() * (1 + axisMargin)),
                    ylim = (intensity[intensity != 0].min()
                                * (1 - axisMargin),
                            intensity.max()
                                * (1 + axisMargin)
                    ),
                    xscale = 'log', yscale = 'log',
                    xlabel = 'q, 1/m', ylabel = 'intensity, 1/(m sr)'
                )
        qAxis = setAxis(qAxis)
        errorbar(q, intensity, intError, zorder = 2, fmt = 'k.',
                 ecolor = 'k', elinewidth = 2, capsize = 4, ms = 5,
                 label = 'Measured intensity', lw = 2,
                 solid_capstyle = 'round', solid_joinstyle = 'miter')
        grid(lw = 2, color = 'black', alpha = .5, dashes = [1, 6],
             dash_capstyle = 'round', zorder = -1)
        aq = numpy.sort(result['fitQ'])
        aI = result['fitIntensityMean'][0, 
                numpy.argsort(result['fitQ'])]
        plot(aq, aI, 'r-', lw = 3, label = 'MC Fit intensity', zorder = 4)
        try:
            plot(aq, numpy.mean(result['scalingFactors'][1, :]) + 0*aq,
                 'g-', linewidth = 3,
                 label = 'MC Background level:\n\t ({0:03.3g})'
                         .format(numpy.mean(result['scalingFactors'][1, :])),
                 zorder = 3)
        except:
            pass
        #for some reason, the axis settings are not set using the above
        #call, and need repeating:
        gca().set_xlim((q.min() * (1 - axisMargin),
                q.max() * (1 + axisMargin)))
        gca().set_ylim((intensity[intensity != 0].min() * (1 - axisMargin),
                            intensity.max() * (1 + axisMargin) ))
        legend(loc = 1, fancybox = True, prop = textfont)
    title('Measured vs. Fitted intensity',
          fontproperties = textfont, size = 'x-large')

    sizeAxis = list()
    # plot histograms
    for parami in range(nHists):
        # get data:
        res = allRes[parami]
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
        plotPar = params[parameterId[parami]]
        plotTitle = plotPar.displayName()
        xLabel = '{}, {}'.format(plotPar.name(), plotPar.suffix())

        # prep axes
        if params[parameterId[parami]].histogram().scaleX == 'log':

            xLim = (histXLowerEdge.min() * (1 - axisMargin), 
                    histXLowerEdge.max() * (1 + axisMargin))
            xScale = 'log'
        else:
            xLim = (histXLowerEdge.min() - (1 - axisMargin)
                    * histXLowerEdge.min(), 
                    histXLowerEdge.max() * (1 + axisMargin))
            xScale = 'linear'

        yLim = (0, volHistYMean.max() * (1 + axisMargin) )
        sizeAxis.append(fig.add_subplot(
                            1, (nHists + 1), parami + 2,
                            axisbg = (.95, .95, .95),
                            xlim = xLim,
                            ylim = yLim,
                            xlabel = xLabel,
                            xscale = xScale,
                            ylabel = '[Rel.] Volume Fraction'))

        sizeAxis[parami] = setAxis(sizeAxis[parami])
        # fill axes
        bar(histXLowerEdge[0:-1], volHistYMean, 
                width = histXWidth, color = 'orange',
                edgecolor = 'black', linewidth = 1, zorder = 2,
                label = 'MC size histogram')
        plot(histXMean, volHistMinReq, 'ro', 
                ms = 5, markeredgecolor = 'r',
                label = 'Minimum visibility limit', zorder = 3)
        errorbar(histXMean, volHistYMean, volHistYStd,
                zorder = 4, fmt = 'k.', ecolor = 'k',
                elinewidth = 2, capsize = 4, ms = 0, lw = 2,
                solid_capstyle = 'round', solid_joinstyle = 'miter')
        legend(loc = 1, fancybox = True, prop = textfont)
        title(plotTitle, fontproperties = textfont,
              size = 'x-large')
        # reapply limits in x
        xlim((histXLowerEdge.min() * (1 - axisMargin),
              histXLowerEdge.max() * (1 + axisMargin)))

    fig.subplots_adjust(left = 0.1, bottom = 0.11,
                        right = 0.96, top = 0.95,
                        wspace = 0.23, hspace = 0.13)
    # trigger plot window popup
    show()

# vim: set ts=4 sts=4 sw=4 tw=0:
