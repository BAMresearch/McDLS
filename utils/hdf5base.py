# -*- coding: utf-8 -*-
# utils/hdf5base.py

""" 
Helper functions for HDF5 functionality
"""

import logging
import h5py

def h5w(wloc, field, hDat, hType = "dataset"):
    """ 
        writes dataset *hDat* to HDF5 location *wloc*, deleting if exists 
        htype can be "dataset" or "attribute"
    """
    # remove old field, only removes link, does not reclaim!
    # ideally, new h5py file should be generated on end:
    # http://stackoverflow.com/questions/11194927/deleting-information-from-an-hdf5-file
    if "dataset" in hType:
        if field in wloc:
            del wloc[field]
        try:
            wloc.create_dataset(field, data = hDat, compression = "gzip")
        except TypeError:
            # a TypeError is raised for non-chunkable data (such as string)
            wloc.create_dataset(field, data = hDat)
        except:
            raise
    else:
        wloc.attrs[field] = hDat

# vim: set ts=4 sts=4 sw=4 tw=0:
