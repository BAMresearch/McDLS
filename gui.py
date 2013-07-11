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

import getopt
from gui import eventLoop

def cmdName(argv):
    if not argv or len(argv) <= 0:
        return ""
    else:
        return os.path.basename(argv[0])

def printUsage(argv, msg=""):
    if len(msg):
        msg += "\n"
    msg += "USAGE: " + cmdName(argv) + " <A> <B> <intensity-file>"
    print >> sys.stdout, msg
    return -1

def main(argv=None):

    if argv is None:
        argv = sys.argv

    opts = []
    args = []
    # extract valid options
    try:
        opts, args = getopt.getopt(argv[1:], "h", ["help"])
    except getopt.error:
        return printUsage(argv, "For a list of valid options, "
                                "use '-h' for help.")
    # evaluate valid options
    if (unicode("-h"), "") in opts or \
       (unicode("--help"), "") in opts: 
        return printUsage(argv)

    return eventLoop()

if __name__ == "__main__":
    main()

# vim: set ts=4 sts=4 sw=4 tw=0:
