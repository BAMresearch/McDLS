#!/usr/bin/env python

import multiprocessing
from .main import main

if __name__ == "__main__":
    multiprocessing.freeze_support()
    # trying to force "spawn" under mac os x
    # multiprocessing.set_start_method('spawn') # not working on Enthought?
    # possbily not existing in python 2.7
    main()

# vim: set ts=4 sts=4 sw=4 tw=0:
