# -*- coding: utf-8 -*-
# mcsas/mcsasparameters.py
# Find the reST syntax at http://sphinx-doc.org/rest.html

import os
import inspect
import logging                                                                 
import json
from utils import units, isString
from utils.propertynames import PropertyNames
from utils.parameter import Parameter
from main import makeAbsolutePath
logging.basicConfig(level = logging.INFO)                                      

class McSASParameters(PropertyNames):
    """
    Defines the static parameters used for the fitting procedure:
        - *model*: an instance of McSAS.model defining the fitting model
        - *contribParamBounds*: Bounds of the active (fitting) parameter
        - *qBounds*: limits in *q* between which the fit is applied
        - *psiBounds*: limits in azimuthal angle (2D patterns only) to which
            the fit is applied
        - *qMagnitude*: indicates the multiplier to scale q to m^-1
        - *iMagnitude*: indicates the multiplier to scale I to (m sr)^-1
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
            datapoints with intensity values < 0
        - *maskZeroInt*: similar to above, but for intensity values = 0

    Most of them should be moved to McSAS as dynamic parameters of type
    *Parameter* which allows them to be configurable in the GUI.
    Some, esp. *histogramBins* and *histogramXScale* should be moved to
    a custom *FitParameter* class along with the *active* flag. (WIP)

    Extension in progress: the class can take keyword-value pairs to 
    overwrite (some of) the parameter values. Additionally, a default
    configuration file (json-style) can be provided using kwarg:
    ``paramDefFile = 'path/to/file'``
    A (limited) set of custom parameters can be supplied in another
    configuration file that overwrite the defaults:
    ``paramFile = 'path/to/file'``

    The order is:

        - default file, whose values are superseded by
        - custom file, whose values are superseded by
        - keyword-value pairs, which only set parameter *values*, or
        - keyword-dict pairs, which sets parameter attributes as defined in
          the supplied dictionary
    """
    # set old-style defaults
    model = None
    contribParamBounds = ()

    # new defaults for loading parameters
    parameters = list()
    paramDefFile = os.path.join("mcsas", "mcsasparameters.json")

    def loadParameters(self, filename):
        # load parameter definitions from file and add to list:
        parDict = dict()
        filename = makeAbsolutePath(filename)
        if not os.path.exists(filename):
            # trying one level up
            dirname, basename = os.path.split(filename)
            dirname = os.path.dirname(dirname)
            filename = os.path.join(dirname, basename)
        try:
            with open(filename, 'r') as jfile:
                logging.info("Loading parameters from file: '{}'"
                             .format(filename))
                parDict = json.load(jfile)
        except IOError:
            logging.error("Could not load default parameter file '{fn}'!"
                          .format(fn = filename))
        for pkey in parDict.keys():
            default = parDict[pkey].pop('default')
            unitClass = parDict[pkey].pop('unitClass', None)
            displayUnit = parDict[pkey].pop('displayUnit', None)
            unitInstance = self.pickUnit(unitClass, displayUnit)

            self.parameters.append(
                    Parameter(pkey, default, unit = unitInstance,
                        **parDict[pkey])
                    )
            logging.debug('Parameter {} ingested'.format(pkey))

    def __init__(self, paramDefFile = None):
        """initialise the defaults and populate the database with values
        where appropriate
        default parameter file can be provided using kwarg:
        paramDefFile = 'path/to/file' relative to application root dir
        """
        # instantiate defaults:
        if not isString(paramDefFile) or not os.path.exists(paramDefFile):
            paramDefFile = self.paramDefFile
        self.loadParameters(paramDefFile)

    def pickUnit(self, unitClass = None, displayUnit = None):
        """ returns a unit object instance of the right class and displayUnit """
        #first make a dictionary of name-object pairs
        if (unitClass is None) or (displayUnit is None):
            logging.warning('Unit not properly set for parameter.')
            return units.NoUnit()

        unitDict = {}
        for name, obj in inspect.getmembers(units):
            if inspect.isclass(obj):
                if issubclass(obj, units.Unit):
                    unitDict.update({name: obj})
        if unitClass in unitDict:
            unitObj = unitDict[unitClass]
            return unitObj(displayUnit)
        else:
            logging.warning('Unit class "{}" not found! Mind the case sensitivity!'
                    .format(unitClass))
            return units.NoUnit()

# vim: set ts=4 sts=4 sw=4 tw=0:
