# -*- coding: utf-8 -*-
# utils/binning.py

"""
 - :py:func:`binningArray`:
   Can be used to do n-by-n pixel binning of 2D detector
   images. The returned uncertainty is the larger of either the binned
   uncertainty or the sample standard deviation in the bin.
 - :py:func:`binning1d`:
   bins the data and propagates errors, or calculates errors
   if not initially provided
 - :py:func:`binningWeighted1d`:
   Weighted binning, where the intensities of a
   pixel are divided between the two neighbouring bins depending on the
   distances to the centres. If error provided is empty, the standard
   deviation of the intensities in the bins are computed.
"""
from __future__ import print_function

from builtins import range
from numpy import (zeros, mean, sqrt, std, reshape, size, linspace,
                   argsort, ones, array, sort, diff)

def binningArray(q, psi, intensity, error, s = 2):
    """This function applies a simple s-by-s binning routine on images.
    It calculates new error based on old error superseded by standard
    deviation in a bin.
    """
    def isodd(x):
        #checks if a value is even or odd
        return bool(x & 1)

    ddi = {'q': q, 'psi': psi, 'intensity': intensity, 'error':error}
    
    sq = q.shape
    if isodd(sq[0]):
        # trim edge
        for it in list(ddi.keys()):
            ddi[it] = ddi[it][1:, :]
    if isodd(sq[1]):
        # trim edge
        for it in list(ddi.keys()):
            ddi[it] = ddi[it][:, 1:]
    # now we can do n-by-n binning
    sq = q.shape
    qo = zeros((sq[0]/s, sq[1]/s))
    ddo = {'q': qo.copy(), 'psi': qo.copy(),
           'intensity': qo.copy(), 'error': qo.copy()}
    for it in 'q', 'psi', 'intensity':
        for ri in range(sq[0]/s):
            for ci in range(sq[1]/s):
                ddo[it][ri, ci] = mean(
                        ddi[it][s*ri:(s*ri+s), s*ci:(s*ci+s)])
        
    for ri in range(sq[0]/s):
        for ci in range(sq[1]/s):
            meanE = sqrt(sum((
                        ddi['error'][s*ri:(s*ri+s), s*ci:(s*ci+s)])**2
                    ))/s**2
            # sample standard deviation
            stdI = std(ddi['intensity'][s*ri:(s*ri+s), s*ci:(s*ci+s)])
            # stdI=0
            ddo['error'][ri, ci] = max((meanE, stdI))
    return ddo

def binning1d(q, intensity, error = None, numBins = 200, stats = 'std'):
    """An unweighted binning routine.
    The intensities are sorted across bins of equal size. If provided error
    is empty, the standard deviation of the intensities in the bins are
    computed.
    """
    
    # Let's make sure the input is consistent
    if size(q) != size(intensity):
        print("Incompatible sizes of q and intensity")
        return
    elif error is not None and size(error) != size(intensity):
        print("Size of error is not identical to q and intensity")
        return

    #flatten q, intensity and error
    q = reshape(q, size(q), 0)
    intensity = reshape(intensity, size(intensity), 0)
    if error is None:
        error = []
    error = reshape(error, size(error), 0)

    # define the bin edges and centres, and find out the stepsize while
    # we're at it. Probably, there is no need for knowing the edges...
    qbinEdges = linspace(q.min(), q.max(), numBins + 1)
    stepsize = qbinEdges[1] - qbinEdges[0]
    qbinCenters = linspace(q.min() + 0.5*stepsize,
                           q.max() - 0.5*stepsize, numBins)
    
    # sort q, let intensity and error follow sort
    sortInd = argsort(q, axis = None)
    q = q[sortInd]
    intensity = intensity[sortInd]
    ibin = zeros(numBins)
    sdbin = zeros(numBins)    
    sebin = zeros(numBins)    
    if size(error) != 0:
        error = error[sortInd]

    # now we can fill the bins
    for bini in range(numBins):
        # limit ourselves to only the bits we're interested in:
        limMask = ((q  > (qbinCenters[bini] - stepsize)) &
                          (q <= (qbinCenters[bini] + stepsize)))
        iToBin = intensity[limMask]
        if size(error) != 0:
            eToBin = sum(error[limMask])

        # find out the weighting factors for each (q, intensity, error)-pair
        # in the array  
        weightFactor = ones(size(iToBin))

        # sum the intensities in one bin and normalize by number of pixels
        ibin[bini] = sum(iToBin)/sum(weightFactor)

        # now we deal with the Errors:
        if (size(error) != 0):
            # if we have errors supplied from outside
            # standard error calculation:
            sebin[bini] = (sqrt(sum(eToBin**2 * weightFactor))/
                            sum(weightFactor))
            if stats == 'auto':
                # according to the definition of sample-standard deviation
                sdbin[bini] = sqrt(sum((
                                iToBin - ibin[bini])**2 * weightFactor
                                   )/(sum(weightFactor) - 1))
                # maximum between standard error and Poisson statistics
                sebin[bini] = array([
                            sebin[bini],
                            sdbin[bini] / sqrt(sum(weightFactor))]).max()
        else:           
            # calculate the standard deviation of the intensity in the bin
            # both for samples with supplied error as well as for those where
            # the error is supposed to be calculated
            # according to the definition of sample-standard deviation
            sdbin[bini] = sqrt(sum(
                                (iToBin-ibin[bini])**2 * weightFactor
                                )/(sum(weightFactor) - 1))
            # calculate standard error by dividing the standard error by the
            # square root of the number of measurements
            sebin[bini] = sdbin[bini]/sqrt(sum(weightFactor))

    return qbinCenters, ibin, sebin

def binningWeighted1d(q, intensity, error = None,
                      numBins = 200, stats = 'se'):
    """Implementation of the binning routine written in Matlab.
    The intensities are divided across the q-range in bins of equal size.
    The intensities of a pixel are divided between the two neighbouring bins
    depending on the distances to the centres. If error provided is empty,
    the standard deviation of the intensities in the bins are computed.

    **Usage**::

        qbin, ibin, ebin = binning_weighted_1d(q, intensity, error = [],
                                               numBins = 200, stats = 'se')

    **Optional input arguments**:

    *numBins*: integer indicating the number of bins to divide the intensity
        over. Alternatively, this can be an array of equidistant bin centres.
        If you go this route, depending on the range, not all intensity may be
        counted.
    *stats*: Can be set to 'auto'. This takes the maximum error between
        supplied Poisson statistics error-based errors or the standard error.

    Written by Brian R. Pauw, 2011, released under BSD open source license.
    """
    # let's make sure the input is consistent
    if size(q) != size(intensity):
        print ("Incompatible lengths of q and intensity, "
               "q must be of the same number of elements as intensity")
        return
    elif error is not None and (size(error) != size(intensity)):
        print("Size of error is not identical to q and intensity")
        return
    stats = stats.lower()
    if (stats != 'std' and stats != 'poisson' and
        stats != 'se' and stats != 'auto'):
        print("Statistics can only be set to 'se' (default), or 'auto'.\n")
        print ("Only use 'auto' for photon-counting detectors, selects "
               "largest error between se and Poisson.\n")
        print ("If errors are supplied, standard errors are not calculated "
               "except in the case of 'auto'")
        return
    if (size(numBins) == 1):
        if (numBins < 1):
            print ("number of bins, numBins, is smaller than one. "
                   "intensity need at least one bin to fill")
            return
    if size(numBins) > 1:
        print ("numBins is larger than one value. Assuming that an "
               "equidistant list of bin centres has been supplied")

    # flatten q, intensity and error
    q = reshape(q, size(q), 0)
    intensity = reshape(intensity, size(intensity), 0)
    if error is None:
        error = []
    error = reshape(error, size(error), 0)
    
    if size(numBins) == 1:
        # define the bin edges and centres, and find out the stepsize while
        # we're at it. Probably, there is no need for knowing the edges...
        dummy, stepsize = linspace(q.min(), q.max(),
                                   numBins + 1, retstep = True)
        qBinCenters = linspace(q.min() + 0.5*stepsize,
                               q.max() - 0.5*stepsize, numBins)
    else:
        if (q.min() > numBins.max() or
            q.max() < numBins.min()):
            print ("Bin centres supplied do not overlap with the q-range, "
                   "cannot continue")
            return
        qBinCenters = sort(numBins)
        stepsize = mean(diff(qBinCenters))
        numBins = size(qBinCenters)

    # initialize output matrices
    ibin = zeros(numBins)
    sdbin = zeros(numBins)    
    sebin = zeros(numBins)    

    # now we can fill the bins
    for bini in range(numBins):
        # limit ourselves to only the bits we're interested in:
        limMask = (q  > (qBinCenters[bini] - stepsize) and
                          q <= (qBinCenters[bini] + stepsize))
        qToBin = q[limMask]
        iToBin = intensity[limMask]
        if (size(error) != 0):
            eToBin = error[limMask]

        # find out the weighting factors for each (q, intensity, error)-pair
        # in the array
        qDist = abs(qToBin - qBinCenters[bini])
        weightFactor = (1 - qDist/stepsize)

        # sum the intensities in one bin
        ibin[bini] = sum(iToBin * weightFactor)/sum(weightFactor)

        # now we deal with the Errors:
        if (size(error) != 0): # if we have errors supplied from outside
            # standard error calculation:
            sebin[bini] = (sqrt(sum(eToBin**2 * weightFactor))
                            / sum(weightFactor))
            if stats == 'auto':
                # according to the definition of sample-standard deviation
                sdbin[bini] = sqrt(sum(
                    (iToBin - ibin[bini])**2 * weightFactor
                                    )/(sum(weightFactor) - 1))
                # maximum between standard error and Poisson statistics
                sebin[bini] = array([
                                sebin[bini],
                                sdbin[bini] / sqrt(sum(weightFactor))]).max()
        else:           
            # calculate the standard deviation of the intensity in the bin
            # both for samples with supplied error as well as for those where
            # the error is supposed to be calculated
            # according to the definition of sample-standard deviation
            sdbin[bini] = sqrt(sum(
                        (iToBin - ibin[bini])**2 * weightFactor
                                )/(sum(weightFactor)-1))
            # calculate standard error by dividing the standard error by the
            # square root of the number of measurements
            sebin[bini] = sdbin[bini]/sqrt(sum(weightFactor))
    return qBinCenters, ibin, sebin
 
# vim: set ts=4 sts=4 sw=4 tw=0:
