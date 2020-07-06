# -*- coding: utf-8 -*-
# utils/loadstore.py

"""
 - :py:func:`pickleLoad`:
   Reads in pickled data from a file (by filename)
 - :py:func:`pickleStore`:
   write a block or dictionary to a file (by filename)
"""

import pickle
from . import mcopen

def pickleLoad(filename):
    """Loads data from a pickle file"""
    fh = mcopen(filename)
    output = pickle.load(fh)
    fh.close()
    return output

def pickleStore(filename, somedata):
    """Writes python object to a file."""
    fh = mcopen(filename, 'w')
    pickle.dump(somedata, fh)
    fh.close()
    return

# vim: set ts=4 sts=4 sw=4 tw=0:
