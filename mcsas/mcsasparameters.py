# -*- coding: utf-8 -*-
# mcsas/mcsasparameters.py
# Find the reST syntax at http://sphinx-doc.org/rest.html

from utils.propertynames import PropertyNames
#more flexible parameter definitions:
from mcsasdefaultcfg import cInfo
import os
import inspect
import logging                                                                 
import json
from utils.parameter import Parameter
logging.basicConfig(level = logging.INFO)                                      
                                                     

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
        - *doPlot*: sets whether to automatically plot or not

       Most of them should be moved to McSAS as dynamic parameters of type
       *Parameter* which allows them to be configurable in the GUI.
       Some, esp. *histogramBins* and *histogramXScale* should be moved to
       a custom *FitParameter* class along with the *active* flag. (WIP)
       
       #######################
       Extension in progress: the class can take keyword-value pairs to 
       overwrite (some of) the parameter values. Additionally, a default
       configuration file (json-style) can be provided using kwarg:
       paramDefFile = 'path/to/file'
       A (limited) set of custom parameters can be supplied in another i
       configuration file that overwrite the defaults: 
       paramFile = 'path/to/file'
       The order is:
       - default file, whose values are superseded by
       - custom file, whose values are superseded by
       - keyword-value pairs, which only set parameter *values*, or
       - keyword-dict pairs, which sets parameter attributes as defined in the
            supplied dictionary

    """
    #set old-style defaults
    model = None
    contribParamBounds = ()
    priors = () # of shape Rrep, to be used as initial guess for
                # analyse(). It will pass on a Prior to MCFit.
    prior = ()  # of shape Rset, to be used as initial guess for
                # MCFit function

    #superseded by new style
    #qBounds = None
    #psiBounds = None
    #deltaRhoSquared = 1.0
    #histogramBins = 50
    #maxRetries = 5
    #maskNegativeInt = False
    #maskZeroInt = False
    #doPlot = False
    #startFromMinimum = False
    #histogramWeighting = 'volume' # can be "volume" or "number"
    #histogramXScale = 'log'
    
    #new defaults for loading parameters
    parameters = list()
    paramDefFile = "McSASParameters.json"

    def loadParameters(self, filename):
        if not os.path.exists(filename):
            logging.error('no default parameter file found!')
            return

        #load parameter definitions from file and add to list:
        with open(filename, 'r') as jfile:
            logging.info('loading parameters from file: {}'.format(filename))
            parDict=json.load(jfile)
        for pkey in parDict.keys():
            default = parDict[pkey].pop('default')
            self.parameters.append(
                    Parameter(pkey, default,
                        **parDict[pkey])
                    )
            logging.debug('Parameter {} ingested'.format(pkey))

    def __init__(self,**kwargs):
        """initialise the defaults and populate the database with values
        where appropriate
        default parameter file can be provided using kwarg:
        paramDefFile = 'path/to/file'
        McSASParameters.json should be in the same directory as this function
        """
        #new style, to gradually replace old style, instantiate defaults:
        fname = self.paramDefFile
        if not(os.path.exists(fname)):
            #try one more:
            #determine the directory in which this module resides
            #determine the directory in which McSASParameters is located:
            #settings should be in the same directory:
            fdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
            fname = os.path.join(fdir, "McSASParameters.json")

        self.loadParameters(fname)

# vim: set ts=4 sts=4 sw=4 tw=0:
