#!/usr/bin/env python

import sys
import os.path

def isFrozen():
    return hasattr(sys, "frozen")

def getScriptPath():
    thisFile = sys.executable
    if not isFrozen():
        thisFile = os.path.join(os.getcwd(), __file__)
    thisFile = os.path.abspath(thisFile)
    return os.path.split(thisFile)

# get script/executable location, add to module search path
scriptPath, scriptFilename = getScriptPath()
programDir = scriptPath
if not isFrozen():
    programDir = os.path.abspath(os.path.join(scriptPath))
    sys.path.append(programDir)

import argparse
import logging
from cutesnake.log import log, replaceStdOutErr

def main(argv = None):
    parser = argparse.ArgumentParser(description = "Monte Carlo SAS analysis")
    parser.add_argument("-t", "--text", action = "store_true",
                        help = "Run in text mode without graphical "
                               "user interface")
    parser.add_argument("-l", "--nolog", action = "store_true",
                        help = "Disable progress output during fit, "
                               "it's written to file in any case.")
    parser.add_argument('fnames', nargs = '*', metavar = 'FILENAME',
                        action = "store",
                        help = "One or more data files to analyse")
    # TODO: add info about output files to be created ...
    args = parser.parse_args()
    from gui.calc import SASData
    SASData.nolog = args.nolog # forwarding logging setting, quick fix for now

    # initiate logging (to console stderr for now)
    replaceStdOutErr() # replace all text output with our sinks

    if not args.text:
        from gui.mainwindow import eventLoop
        # run graphical user interface
        return eventLoop()
    else:
        try:
            from gui import calc
            calc(args.fnames)
        except StandardError, e:
            # show detailed error traceback if there was one
            import traceback
            logging.error(traceback.format_exc())

if __name__ == "__main__":
    main()

# vim: set ts=4 sts=4 sw=4 tw=0:
