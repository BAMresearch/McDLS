#!/usr/bin/env python                                                          
#coding: utf8    
"""
carrier for settings  and info used for a McSAS run
The instance can be imported, read and modified in a module after
calling "import McSASCfg", f.ex. in the module's "__init__"

The parameters are then available in 
McSASCfg.parameters.parameters.[parname]
which can be mapped to more convenient names.
Alternatively, McSASCfg.parameter.getPar([parname]) returns the
handle to a parameter.
The parameter definitions are read from the McSASParameters.json
dictionary, and new parameters can be easily defined through there.
"""

__author__ = "Brian R. Pauw"
__contact__ = "brian@stack.nl"
__license__ = "GPLv3+"
__copyright__ = "National Institute of Materials Science, Tsukuba, Japan"
__date__ = "2013-12-21"
__status__ = "alpha"
version = "0.0.1"

from cutesnake.algorithm import Parameter
import numpy as np
import logging
import os
from McSASDefaultCfg import cInfo

# instantiate a parameters class filled with defaults
parameters = cInfo()


# this could be a useful function:
def saveSettings(self,filename):
    """save the current parameters to a file"""
    pass

