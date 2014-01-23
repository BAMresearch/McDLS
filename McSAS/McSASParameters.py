# -*- coding: utf-8 -*-
# McSAS.py
# Find the reST syntax at http://sphinx-doc.org/rest.html

from utils.propertynames import PropertyNames

class McSASParameters(PropertyNames):
    """Defines the static parameters used for the fitting procedure:
        - *model*: an instance of McSAS.model defining the fitting model
        - *contribParamBounds*: Bounds of the active (fitting) parameter
        - *qBounds*: limits in *q* between which the fit is applied
        - *psiBounds*: limits in azimuthal angle (2D patterns only) to which
            the fit is applied
        - *priors*: a set of contribution parameters for a number of McSAS 
            repetitions. These can be used to resume or recalculate an older
            result
        - *prior*: same as *priors*, but only for a single repetition. Can be
            used as an initial guess. WARNING: intial guesses may skew result
        - *histogramBins*: number of bins to use for size distribution
            determination. Does not affect fit
        - *histogramXScale*: can be "log" or "linear", sets the horizontal
            axis scaling. Does not affect fit
        - *histogramWeighting*: can be "volume"(recommended) or "number". 
            Determines whether to plot the volume-weighted size distribution
            (which is closest to what is measured in a scattering measurement)
            or a number-weighted size distribution.
        - *deltaRhoSquared* the value (in m^{-4}) of the squared electron 
            density contrast. Typically on the order of 10^{25} to 10^{30} for
            X-ray scattering measurements. To get a correct volume fraction, 
            the intensity must be in reciprocal meters, as must *q* and the
            object size parameters must be adjusted accordingly.
        - *startFromMinimum*: may be depreciated: starts with an initial guess
            consisting of minimum size values rather than random. For testing
            purposes only.
        - *maxRetries*: if a single optimisation fails, it will be retried 
            this integer of times.
        - *maskNegativeInt*: may be depreciated due to overlap with similar 
            functionality in the McSAS.data molule. Setting this will ignore 
            datapoints with intensity values <0
        - *maskZeroInt*: similar to above, but for intensity values =0
        - *lowMemoryFootprint*: can be used for models which require a lot 
            of memory to calculate. Slows down fitting by about 30%, but 
            reduces the memory requirements drastically as only a single 
            contribution is calculated at a time
        - *doPlot*: sets whether to automatically plot or not

       Most of them should be moved to McSAS as dynamic parameters of type
       *Parameter* which allows them to be configurable in the GUI.
       Some, esp. *histogramBins* and *histogramXScale* should be moved to
       a custom *FitParameter* class along with the *active* flag. (WIP)

    """
    model = None
    contribParamBounds = ()
    qBounds = None
    psiBounds = None
    priors = () # of shape Rrep, to be used as initial guess for
                # analyse(). It will pass on a Prior to MCFit.
    prior = ()  # of shape Rset, to be used as initial guess for
                # MCFit function
    histogramBins = 50
    histogramXScale = 'log'
    histogramWeighting = 'volume' # can be "volume" or "number"
    deltaRhoSquared = 1.0
    startFromMinimum = False
    maxRetries = 5
    maskNegativeInt = False
    maskZeroInt = False
    lowMemoryFootprint = False
    doPlot = False

# vim: set ts=4 sts=4 sw=4 tw=0:
