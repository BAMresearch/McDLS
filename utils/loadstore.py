# -*- coding: utf-8 -*-
# utils/loadstore.py

"""
 - :py:func:`pickleLoad`:
   Reads in pickled data from a file (by filename)
 - :py:func:`pickleStore`:
   write a block or dictionary to a file (by filename)
"""

import pickle

def pickleLoad(filename):
    """*\*args* can be 1-4, indicates number of output variables.
    If it is even possible to extract more from pickle."""
    fh = open(filename)
    output = pickle.load(fh)
    fh.close()
    return output

def pickleStore(filename, somedata):
    """Writes DBlock to a file."""
    fh = open(filename, 'w')
    pickle.dump(somedata, fh)
    fh.close()
    return

# vim: set ts=4 sts=4 sw=4 tw=0:
