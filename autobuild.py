#!/usr/bin/env python
# automated build script to be run by auto start on each system
# supposed to be compatible with multiple platforms,
# currently supported and tested: Win7

import sys
import subprocess
import os.path
import logging

logging.basicConfig(level = logging.DEBUG)

def waitForUser():
    print "hit <enter> to continue"
    sys.stdout.flush()
    input = sys.stdin.readline()

WORKDIR = os.path.expanduser('~')
logging.info("Work dir: '{WORKDIR}'".format(**locals()))
REPOURL = "https://bitbucket.org/pkwasniew/mcsas.git"
REPOURL = os.path.abspath("z:\mcsas")
DIRNAME = os.path.basename(os.path.splitext(REPOURL)[0])
logging.info("Repository location: '{REPOURL}'".format(**locals()))

def git(*args, **kwargs):
    cmd = ('git',) + args
    try:
        if not os.path.exists(WORKDIR):
            raise OSError("The work dir '{WORKDIR}' does not exist!"
                          .format(**globals()))
        p = subprocess.Popen(cmd, cwd = WORKDIR,
                             stdout = subprocess.PIPE,
                             stderr = subprocess.PIPE)
    except OSError, e:
        logging.error("GIT was not found! "
                      "Please ensure it is installed and in $PATH!")
        logging.error(e)
        cleanup(1)
    out, err = p.communicate()
    out = out.strip()
    err = err.strip()
    cmd = " ".join(cmd)
    #print >>sys.__stderr__, cmd, p.returncode, out, err
    if p.returncode != kwargs.get("expectedReturnCode", 0):
        logging.error(err)
        cleanup(1)
    return cmd, "\n".join([s for s in out, err if len(s) > 0])

def cleanup(exitCode):
    waitForUser()
    sys.exit(exitCode)

def testForGit():
    git(expectedReturnCode = 1)

def clone():
    if os.path.exists(os.path.join(WORKDIR, DIRNAME)):
        return
    cmd, out = git("clone", REPOURL)
    logging.info("\n".join((cmd, out)))

def getLatestBranch():
    cmd, out = git("for-each-ref", "--sort=-committerdate", "--count=1",
                   "--format=%(refname:short)", "refs")
    return os.path.basename(out.strip())

def getDateTime():
    cmd, out = git("log", "-n1", "--pretty=%ci")
    fields = out.strip().split()
    return "_".join(fields[:-1]).replace(":", "-")

def getHash():
    cmd, out = git("log", "-n1", "--pretty=%h")
    return out.strip()

def checkout(refname):
    logging.info("\n".join(git("checkout", refname)))

def freeze(*args):
    cmd = ("python", "cxfreeze.py", "build_exe") + args
    p = subprocess.Popen(cmd, cwd = WORKDIR,
                         stdout = subprocess.PIPE,
                         stderr = subprocess.PIPE)
    out, err = p.communicate()
    out = out.strip()
    err = err.strip()
    cmd = " ".join(cmd)
    print out
    print err
    print cmd

testForGit()
clone()
WORKDIR = os.path.join(WORKDIR, DIRNAME)
branch = getLatestBranch()
checkout(branch)
datetime = getDateTime()
hash = getHash()
basename = "{DIRNAME} {datetime} {branch} {hash}".format(**locals())
logging.info(basename)
freeze(datetime, branch, hash)

waitForUser()

