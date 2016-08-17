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

# from utils.devtools import DBG

class FindModels(object):
    """ 
    Finds all methods of type ScatteringModel in the subdirectories starting 
    from startPwd. startPwd defaults to the root mcsas pwd + "models".
    returns a list of full paths, and a list of associated model names
    """
    _startPwd = None
    _modelFiles = []
    
    def __init__(self, startPwd = None):
        self._startPwd = startPwd 
        if self._startPwd is None:
            # find root path of main.py, use that: TODO: abspath necessary?
            self._startPwd = os.path.dirname(os.path.abspath(inspect.stack()[-1][0]) )
            self._startPwd = os.path.join(self._startPwd, "models")
            logging.info("Start path for model finding not supplied, using default: {}"
                    .format(self._startPath))
        for filename in self.candidateFiles():
            if isModel(filename):
                self._modelFiles.append(filename)


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

    def isModel(self, modelPath):
        """ 
        Tests whether the functions in a file is a ScatteringModel. 
        Functionality redundant due to modelClassName
        """
        # first test: does it mention "ScatteringModel" in the file? Superseded by modelClassName
        # strFound = False
        # for fline in open(modelPath):
        #     if "SASModel" in fline:
        #         strFound = True
        # if not strFound:
        #     return False
        # second test, does it define a method with the same name as the filename?
        modelClassName = modelPath
        if modelClassName is None:
            return False

        return True

    def modelClassName(modelPath):
        """
        extract the model class name from the module file
        """
        searchPrepend = "class "
        searchAppend = "(sasmodel"
        filebasename = os.path.basename(modelPath)[:-3] # we know it ends with ".py"
        for fline in open(modelPath):
            if searchPrepend + filebasename + searchAppend in fline.lower():
                starti = fline.lower().index(filebase)
                return fline[starti : starti + len(filebasename) + 1]
        return None



    def prioritize(self):
        """ 
        prioritises the lists to put a few primary models on top, like spheres, cylinders
        and ellipsoids. These are predefined.
        """
        pass
    
# vim: set ts=4 sts=4 sw=4 tw=0:
