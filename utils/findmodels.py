# -*- coding: utf-8 -*-
# utils/modelfinder.py

"""
For use in gui/modelwidget.py, to help find valid calculation models
"""

import logging
import inspect
import os
import sys
import imp
from collections import OrderedDict
from models.scatteringmodel import ScatteringModel
from main import makeAbsolutePath
from utils import isList

def _getValidClasses(mod, isTopLevel = True):
    if not inspect.ismodule(mod):
        return []
    try:
        if not isTopLevel and not mod.__isModelGroup__:
            return []
    except AttributeError:
        return []
    validClasses = []
    for attrName in dir(mod):
        attr = getattr(mod, attrName)
        if inspect.ismodule(attr):
            validClasses += _getValidClasses(attr, isTopLevel = False)
            continue
        if not inspect.isclass(attr):
            continue
        try: # ignore stuff defined elsewhere, imported into this module
            if attr.__module__ != mod.__name__:
                continue
        except AttributeError:
            continue
        if not issubclass(attr, ScatteringModel):
            continue # focus on scattering models
        if attr.name() is None or "dummy" in attrName.lower():
            continue # ignore dummy models, e.g. from testing code
        validClasses.append(attr)
    return validClasses

def findFiles(searchPath, extension):
    """
    generator for files ending in .py
    code from: http://stackoverflow.com/questions/2186525
    """
    return [ os.path.join(dirpath, filename)
            for dirpath, dirnames, filenames in os.walk(searchPath)
                for filename in filenames
                    if filename.endswith(extension)
                        and not filename.startswith("__init__")]

def reorder(indata, priorityKeys):
    def getMatchingKeys(inmap, partialKey):
        return [key for key in inmap.keys() if partialKey in key]
    outdata = OrderedDict()
    # copy over the entries with priority keys
    for name in priorityKeys:
        keys = getMatchingKeys(indata, name)
        if not len(keys):
            continue
        key = keys[0]
        outdata[key] = indata.pop(key)
    # finally get the remaining entries
    outdata.update(indata)
    return outdata

class FindModels(object):
    """
    Finds all methods of type ScatteringModel in the subdirectories starting
    from searchPath. searchPath defaults to the root mcsas pwd + "models".
    returns a list of full paths, and a list of associated model names
    """
    _searchPath = None
    _modelFiles = {}
    _priorityModels = ( # ordered list of default models appear first
            "Sphere",
            "LMADenseSphere",
            "EllipsoidsIsotropic",
            "CylindersIsotropic",
            "SphericalCoreShell",
            "EllipsoidalCoreShell",
            "GaussianChain",
            "Kholodenko"
    )
    _rootName = "models"
    _models = None

    @property
    def rootName(self):
        return self._rootName

    def __init__(self, searchPath = None):
        self._models = OrderedDict()
        self._searchPath = searchPath
        if self._searchPath is None:
            # base all relative paths on directory of main.py
            self._searchPath = makeAbsolutePath(self.rootName)
            logging.info("Start path for model finding not supplied, "
                         "using default: '{}'".format(self._searchPath))
        for filename in findFiles(self._searchPath, '.py'):
            validModelClasses = self._getModelClasses(filename)
            if not isList(validModelClasses):
                continue
            for cls in validModelClasses:
                # get full class hierachy name, use that as key
                # there may be multiple models per file/module
                key = '.'.join((cls.__module__, cls.__name__))
                self._models[key] = (cls, filename)
        self._models = reorder(self._models, self._priorityModels)
        logging.debug("Ordered model files: {}".format(self._models))

    def _getTailDirs(self, dirpath):
        """Returns the directory names following the self._searchPath base
        path as reversed list.
        Returns an empty list if the given *dirpath* equals the search path.
        """
        rel = os.path.relpath(dirpath, self._searchPath).strip()
        if not len(rel) or rel == ".":
            return []
        lst = []
        while len(rel):
            rel, tail = os.path.split(rel)
            lst.append(tail)
        return reversed(lst)

    def _getModelClasses(self, filePath):
        """extract the model class name from the module file
        """
        dirpath, filename = os.path.split(filePath)
        moduleName, ext = os.path.splitext(filename)
        if "scatteringmodel" in moduleName:
            # do not attempt to reload ScatteringModel class
            # otherwise issubclass() checks will fail (other obj ID)
            return []
        modName = self.rootName # base name for all models full package name
        # treat each sub dir found as module, import if necessary
        for name in self._getTailDirs(dirpath):
            modName += "." + name
            if modName not in sys.modules:
                sys.modules[modName] = imp.new_module(modName)
        moduleName = modName + "." + moduleName
        dirpath = [ os.path.dirname(self._searchPath) ]
        mod = None
        for modName in moduleName.split('.'):
            try:
                result = imp.find_module(modName, dirpath)
            except ImportError:
                logging.warning("Could not load '{}' from {}! "
                                .format(modName, dirpath)
                                + "__init__.py missing?")
#                raise # for debugging, comment for proper operation
                return []
            if mod is not None:
                modName = '.'.join((mod.__name__, modName))
            # loads the file as standalone module
            subMod = imp.load_module(modName, *result)
            if mod is not None:
                subMod.__package__ = mod.__name__
            if hasattr(subMod, "__path__"):
                dirpath = subMod.__path__
            mod = subMod
        return _getValidClasses(mod)

    def __len__(self):
        return len(self._models)

    def __getitem__(self, key):
        return self._models[key]

    def __iter__(self):
        return iter(self._models.keys())

# vim: set ts=4 sts=4 sw=4 tw=0:
