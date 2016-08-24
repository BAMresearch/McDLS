# -*- coding: utf-8 -*-
# utils/modelfinder.py

"""
For use in gui/modelwidget.py, to help find valid calculation models
"""

import logging
import inspect
import os
import h5py
import fnmatch
from models.scatteringmodel import SASModel
from collections import OrderedDict


# from utils.devtools import DBG

class FindModels(object):
    """ 
    Finds all methods of type ScatteringModel in the subdirectories starting 
    from startPwd. startPwd defaults to the root mcsas pwd + "models".
    returns a list of full paths, and a list of associated model names
    """
    _startPwd = None
    _modelFiles = {}
    _priorityModels = [ # list of standard models in order of appearance
            "Sphere",
            "LMADenseSphere",
            "EllipsoidsIsotropic",
            "CylindersIsotropic",
            "SphericalCoreShell",
            "EllipsoidalCoreShell",
            "GaussianChain",
            "Kholodenko"
            ]
    _orderedModelFiles = OrderedDict()
    
    def __init__(self, startPwd = None):
        self._startPwd = startPwd 
        if self._startPwd is None:
            # find root path of main.py, use that: TODO: abspath necessary?
            self._startPwd = os.path.dirname(os.path.abspath(inspect.stack()[-1][1]) )
            self._startPwd = os.path.join(self._startPwd, "models")
            logging.info("Start path for model finding not supplied, using default: {}"
                    .format(self._startPwd))
        for filename in self.candidateFiles():
            modelClassName = self.modelClassName(filename)
            if modelClassName is not None:
                self._modelFiles.update(
                        { self.modelClassName(filename): filename })
        self.prioritize()
        logging.debug("Ordered model files: {}".format(self._orderedModelFiles))

    def candidateFiles(self):
        """
        generator for files ending in .py
        code from: http://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python
        """
        pattern = '*.py'
        for r, d, f in os.walk(self._startPwd):
            for filename in f:
                if fnmatch.fnmatch(filename, pattern):
                    filename = os.path.join(r, filename)
                    yield filename

    def modelClassName(self, modelPath):
        """
        extract the model class name from the module file
        """
        searchPrepend = "class "
        searchAppend = "(sasmodel"
        filebasename = os.path.basename(modelPath)[:-3] # we know it ends with ".py"
        for fline in open(modelPath):
            if searchPrepend + filebasename + searchAppend in fline.lower():
                starti = fline.lower().index(filebasename)
                return fline[starti : starti + len(filebasename)]
        return None

    def prioritize(self):
        """ 
        prioritises the lists to put a few primary models on top, like spheres, cylinders
        and ellipsoids. These are predefined. Returns an orderedDict of model names and paths
        """
        modelNames = self._modelFiles.keys()
        # first we fill in the models matching the prioritized models
        for oName in self._priorityModels:
            if oName in modelNames:
                modelNames.remove(oName)
                self._orderedModelFiles.update({oName: self._modelFiles[oName]})

        # remainder of names are listed alphabetically. 
        modelNames = sorted(modelNames)
        for mName in modelNames:
            self._orderedModelFiles.update({mName: self._modelFiles[mName]})

    
# vim: set ts=4 sts=4 sw=4 tw=0:
