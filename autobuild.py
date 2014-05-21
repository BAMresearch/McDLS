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

import re
from distutils.version import StrictVersion
import requests
assert(StrictVersion(requests.__version__) > StrictVersion('2.2.0'))
from collections import OrderedDict

class BitbucketClient(object):
    """Minimal Bitbucket that signs in and uploads files.
    from https://bitbucket.org/dahlia/bitbucket-distutils
    version: 0.1.2 (public domain, 2014-05-20)
    """

    _DBG = False
    _contentFn = "last_content.html"
    _credentialsFn = "credentials"

    def __init__(self, repository, username = None, password = None):
        if username is None or password is None:
            username, password = self.getCredentials(self._credentialsFn)
        if not self._DBG:
            self.session = requests.session()
            self.signin(username, password)
        self.repository = repository

    def signin(self, username, password):
        url = 'https://bitbucket.org/account/signin/'
        form = self.session.get(url)
        token = self._find_field(form.content, 'name', 'csrfmiddlewaretoken')
        data = {'username': username, 'password': password,
                'csrfmiddlewaretoken': token}
        login = self.session.post(url, data = data, cookies = form.cookies,
                                  headers = {'Referer': url})
        self.cookies = login.cookies

    def upload(self, filename):
        assert(os.path.exists(filename))
        url = 'https://bitbucket.org/' + self.repository + '/downloads'
        fields = ('acl', 'success_action_redirect', 'AWSAccessKeyId',
                  'Policy', 'Signature', 'Content-Type', 'Content-Disposition',
                  'key',)
        if self._DBG:
            with open(self._contentFn) as fd:
                content = fd.read()
        else:
            form = self.session.get(url, cookies = self.cookies)
            content = form.content
            with open(self._contentFn, 'w') as fd:
                fd.write(content)
        s3_url = self._find_field(content, "id", "upload-file-form",
                                  ftype = "form", fquery = "action")
        data = OrderedDict((f, self._find_field(content, "name", f)) for f in fields)
        basename = os.path.basename(filename)
        data['Content-Type'] = "application/octet-stream"
        data['Content-Disposition'] += basename
        data['key'] += basename
        if self._DBG:
            print data
            return None
        with open(filename, 'rb') as fp:
            files = {'file': (basename, fp)}
            response = self.session.post(s3_url, data = data, files = files)
        if 300 <= response.status_code < 400 and 'location' in response.headers:
            response = self.session.get(response.headers['location'])
            logging.warning("Updated response!")
        logging.info("Response: {code} {hdr}"
                     .format(code = response.status_code,
                             hdr = str(response.headers)))
        assert(200 <= response.status_code < 300)
        return url + '/' + basename

    def _find_field(self, content, key, value,
                    ftype = "input", fquery = "value"):
        """Extracts a given query value from an HTML tag with a given field
        type and key/value identification."""
        pattern = (r'<' + ftype + '\s[^<>]*'
                        + key + '=[\'"]' + re.escape(value)
                        + r'[\'"]\s[^>]*>')
        tag = re.search(pattern, content)
        token = re.search(fquery + r'=[\'"]([^\'"]*)[\'"]',
                          tag.group(0))
        return token.group(1)

    @classmethod
    def getCredentials(cls, filename):
        programDir = os.path.dirname(
                        os.path.abspath(
                            os.path.join(os.getcwd(), __file__)))
        filename = os.path.join(programDir, filename)
        logging.info("Reading {cls} credentials from '{fn}'!"
                     .format(cls = cls.__name__, fn = filename))
        creds = ""
        with open(filename) as fd:
            creds = fd.read()
        return creds.split()

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

# TODO: let freeze output the created file, upload it below
# TODO further:
#   - setup logging into a file, outside with python stderr eventually?
#     -> push it to a special branch ('buildlogs' or similar)
#     -> make it available via https://readthedocs.org/projects/mcsas/
#        (allows some formatting)
#   - same goes for reports of automated tests later
#     -> will be plugged in right after cloning

#bb = BitbucketClient("pkwasniew/mcsas")
#print bb.upload(os.path.abspath("README.md"))
#url = self.upload("/home/ingo/code/mcsas/McSASGui-2014-05-19_19-47-53-restructuring-157d5b5.7z")

waitForUser()

# vim: set ts=4 sts=4 sw=4 tw=0:
