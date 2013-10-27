# -*- coding: utf-8 -*-
# lastpath.py

"""
General utilities without GUI dependencies.
"""

import os
import logging
from cutesnake.utils import isString

def getHomeDir():
    return os.path.expanduser('~')

class LastPath(object):
    """
    Stores a file system path for use in file open dialogs.

    How to test this platform independent?
    >>> from utils import LastPath, getHomeDir
    >>> LastPath.path == getHomeDir()
    True
    >>> LastPath.path = '.'
    >>> LastPath.path == '.'
    True
    """

    _path = getHomeDir()

    @classmethod
    def get(cls):
        return cls._path

    @classmethod
    def set(cls, lastpath):
        """Accepts a directory path or a file path.
        Determines the directory itself in the latter case."""
        if not isString(lastpath):
            lastpath = unicode(lastpath)
        # get path of possible unwritten files (previously selected)
        path = lastpath
        while not os.path.isdir(path) and len(path) > 0:
            path = os.path.dirname(path)
        if len(path) > 0:
            cls._path = path
        else:
            logging.warning("Could not set last path '{0}'!".format(lastpath))

# vim: set ts=4 sts=4 sw=4 tw=0:
