#!/usr/bin/env python                                                          
#coding: utf8    
"""
default for settings  and info used for a McSAS run
used by McSASCfg
"""

__author__ = "Brian R. Pauw"
__contact__ = "brian@stack.nl"
__license__ = "GPLv3+"
__copyright__ = "National Institute of Materials Science, Tsukuba, Japan"
__date__ = "2013-12-21"
__status__ = "alpha"
version = "0.0.1"

from utils.parameter import (Parameter, ParameterFloat,
        ParameterBoolean, ParameterNumerical, ParameterString)
import logging, json
import os, inspect
import numpy as np

class ExtendedEncoder(json.JSONEncoder):
    """JSON encoder extended to deal with Unicode, arrays and cls descriptions"""
    def default(self, obj):
        if isinstance(obj, unicode):
            return str(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif obj == np.array or obj == np.ndarray :
            return "array"
        elif isinstance(obj, type):
            #return a suitable string
            if obj is float:
                return "float"
            elif obj is unicode:
                return "unicode"
            elif obj is str:
                return "str"
            elif obj is int:
                return "int"
            else:
                return repr(obj)

        return json.JSONEncoder.default(self, obj)

class cInfo(object):
    """
    This class contains all the information required to read, verify and write
    configuration parameters files.
    """
    parameters=None
    logging.getLogger('McSAScfg')
    logging.basicConfig(level = logging.DEBUG)
    parameterNames=list()

    def __init__(self,**kwargs):
        """initialise the defaults and populate the database with values
        where appropriate
        default parameter file can be provided using kwarg:
        paramDefFile = 'path/to/file'
        McSASParameters.json should be in the same directory as this function
        """
        fname = kwargs.get("paramDefFile", None)
        if fname is None:
            if os.path.exists("McSASParameters.json"):
                fname = "McSASParameters.json"
            else:
                #try one more:
                #determine the directory in which McSASDefaultsCfg is located:
                #settings should be in the same directory:
                fdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
                fname = os.path.join(fdir, "McSASParameters.json")

        if not os.path.exists(fname):
            logging.error('no default parameter file found!')
            return false

        self.loadParams(fname = fname)

    def loadParams(self, fname = None):
        """
        writes the default definitions and bounds for the configuration
        parameters to self.parameters
        Can also be used to update existing parameters from supplied filename
        """
        with open(fname, 'r') as jfile:
            logging.info('loading parameters from file: {}'.format(fname))
            parDict=json.load(jfile)

        if self.parameters is None:
            #create if it does not exist yet
            self.parameters=lambda: None

        #now we cast this information into the Parameter class:
        for kw in parDict.keys():
            subDict = parDict[kw]
            name = kw
            value = subDict.pop("value", None)
            default = subDict.pop("default", None)
            if value is None and default is not None:
                value = default
            #determine parameter class:
            cls = subDict.pop("cls", None)
            if cls == "int":
                subDict.update(cls = ParameterNumerical)
            elif cls == "float":
                subDict.update(cls = ParameterFloat)
            elif cls == "bool": 
                subDict.update(cls = ParameterBoolean)
            elif cls == "str": 
                subDict.update(cls = ParameterString)
            else:
                logging.warning('parameter type {} for parameter {} not understood from {}'.format(cls, kw, fname ))

            if kw in self.parameterNames:
                #value exists, should be updated with supplied kwargs
                self.set(kw,**subDict)
                logging.info('successfully updated parameter: {}'.format(kw))
            else:
                #logging.info('ingesting parameter {} with value {}'.format(name, value))
                temp = Parameter(name, value, **subDict)
                setattr(self.parameters,kw,temp)
                self.parameterNames.append(kw)
                logging.info('successfully ingested parameter: {}'.format(kw))
        

    def writeConfig(self, fname):
        """
        writes the configuration to a settings file. 
        Required input parameter is the filename to write to.
        """

        #create a dictionary containing one (sub-)dictionary per parameter
        parDict = dict()
        for kw in self.parameterNames:
            subDict = dict()
            par = self.getPar(kw)
            for aName in par.attributeNames():
                subDict[aName] = par.get(aName)

            parDict[kw] = subDict

        #write dictionary to file
        with open(fname,'wb') as jfile:
            json.dump(parDict, jfile, cls=ExtendedEncoder, indent = 4, sort_keys=True)


    def parseConfig(self):
        """
        Runs through the entire settings, raising warnings where necessary
        """
        for pn in self.parameterNames:
            pf = self.getPar(pn)
            if pf.value is None:
                continue #skip this value, has not yet been set
            # rewrite for McSAS:
            #pf.checkSize()
            #pf.clipValue()
            
    def getPar(self,key):
        #returns the handle to the parameter defined by key or returns None
        #if it doesn't exist
        if key in self.parameterNames:
            return getattr(self.parameters,key)
        else:
            logging.warning(
                    'Could not find parameter {}. Define base parameter in parameter settings file first'.format(key)
                    )
            return None

    def setParVal(self,par,value):
        "shortcut method for setting the value of a parameter only"
        parhandle = self.getPar(par)
        parhandle.setValue(value)
    
    def getParVal(self,par):
        "shortcut method for getting the value of a parameter"
        parhandle = self.getPar(par)
        return parhandle.value()

    #def set(self,par,**kwargs):
    #    """
    #    sets one or more parameter attributes
    #    not sure the "set" function works in Parameter as in imp2
    #    TODO: update for McSAS
    #    """
    #    parhandle = self.getPar(par)
    #    
    #    for kw in kwargs:
    #        parhandle.set(kw, kwargs[kw])


