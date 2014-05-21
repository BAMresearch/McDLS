#!/usr/bin/env python

import sys
import os.path

def getScriptPath():
    """Returns the full path to the current script file which calls this
    function. See cutesnake.utils"""
    thisFile = sys.executable # works for frozen app
    try: # __file__ not defined if frozen
        thisFile = os.path.join(os.getcwd(), __file__)
    except NameError: pass
    thisFile = os.path.abspath(thisFile)
    return os.path.split(thisFile)

# get script/executable location, add to module search path
SCRIPT_PATH, SCRIPT_FILENAME = getScriptPath()
if not hasattr(sys, "frozen") and SCRIPT_PATH not in sys.path:
    sys.path.append(SCRIPT_PATH)

def makeAbsolutePath(relpath):
    return os.path.abspath(os.path.join(SCRIPT_PATH, relpath))

import argparse
import logging
from cutesnake.log import replaceStdOutErr
import gui.calc

def main(argv = None):
    parser = argparse.ArgumentParser(description = "Monte Carlo SAS analysis")
    parser.add_argument("-t", "--text", action = "store_true",
                        help = "Run in text mode without graphical "
                               "user interface")
    parser.add_argument("-l", "--nolog", action = "store_true",
                        help = "Disable progress output during fit, "
                               "it's written to file in any case.")
    parser.add_argument("-s", "--start", action = "store_true",
                        help = "Start the calculation immediately.")
    parser.add_argument('fnames', nargs = '*', metavar = 'FILENAME',
                        action = "store",
                        help = "One or more data files to analyse")
    # TODO: add info about output files to be created ...
    args = parser.parse_args()
    # forwarding logging setting, quick fix 
    gui.calc.Calculator.nolog = args.nolog

    # initiate logging (to console stderr for now)
    replaceStdOutErr() # replace all text output with our sinks

    if not args.text:
        from gui.mainwindow import eventLoop
        # run graphical user interface, passing argument parser result
        return eventLoop(args)
    else:
        # TODO: fix command line run
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
