# -*- coding: utf-8 -*-
# utils/hdf5base.py

""" 
Helper functions for HDF5 functionality
"""

import logging
import h5py
from abc import ABCMeta

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

class HDF5Mixin(object):
    __metaclass__ = ABCMeta
    _h5Datasets = []
    _h5Attrs = []
    _h5Callers = []
    _h5LocAdd = ''

    # mixin object (?) to add a write function to other classes
    def __init__(self, **kwargs):
        super(HDF5Mixin, self).__init__(**kwargs)

    # @mixedmethod
    # @classmethod
    def writeHDF(self, filename, loc, item = None):
        """ 
        Writes the vector to an HDF5 output file *filename*, at location *loc*. 
        This location should be e.g. "/mcentry01/[ sas | dls ]data01/x0". 
        """
        loc += self._h5LocAdd
        if item is not None:
            loc += "/" + item
        # process datasets and attribute calls
        logging.debug("Writing at loc: {}".format(loc))
        logging.debug("_h5Datasets: {}, _h5Attrs: {}"
                .format(self._h5Datasets, self._h5Attrs))
        with h5py.File(filename) as h5f:
            wloc = h5f.require_group(loc) # unicode's no problem
            for field in self._h5Datasets:
                hDat = getattr(self, field, None)
                if hDat is not None:
                    logging.debug("*** dataset name: {}, data: {}"
                            .format(field, hDat))
                    h5w(wloc, field, hDat, hType = "dataset")
            for att in self._h5Attrs:
                logging.debug("*** attribute field: {}, data: {}"
                        .format(att, self.name))
                h5w(h5f[loc], att, self.name, hType = "attribute")

        # process custom callers
        for call in self._h5Callers:
            caller = getattr(self, call, None)
            if caller is None:
                # skip it
                continue
            callerWFunc = getattr(caller, "writeHDF", None)
            if callerWFunc is not None:
                logging.debug("Calling HDF5 writer for {}, at loc: {}"
                        .format(call, loc))
                callerWFunc(filename, loc, item = call)
            else:
                logging.warning("item {} does not have writeHDF functionality"
                        .format(call))


