#!/usr/bin/env python2
# prepares an iconset according to https://stackoverflow.com/a/20703594
# run it in directory 'resources/icon' to create an iconset directory
# for use with iconutil on OSX:
# 'iconutil -c icns mcsas.iconset' creates mcsas.icns
# requires png images for each required resolution to exist

import os
import sys
import glob
import shutil
from collections import OrderedDict

def resolutions():
    return 16, 32, 128, 256, 512

if __name__ == "__main__":
    # find all required resolutions
    res = dict.fromkeys(resolutions())
    for r in res.keys():
        res[r] = (r, r*2)
    # find all required files
    path = os.getcwd()
    files = dict()
    for r in sum(res.values(), ()):
        fn = glob.glob(os.path.join(path, "*{}.png".format(r)))
        assert len(fn)
        fn = fn[0]
        assert os.path.isfile(fn)
        files[r] = fn
    # create a new directory for the iconset
    targetDir = os.path.join(path, "mcsas.iconset")
    os.mkdir(targetDir)
    # copy source images to their appropriate places&names
    for r, srcRes in res.iteritems():
        basename = "icon_{0}x{0}".format(r)
        assert len(srcRes) == 2
        srcFn = files[srcRes[0]]
        dstFn = os.path.join(targetDir, basename + ".png")
        shutil.copyfile(srcFn, dstFn)
        srcFn = files[srcRes[1]]
        dstFn = os.path.join(targetDir, basename + "@2x.png")
        shutil.copyfile(srcFn, dstFn)

# vim: set ts=4 sts=4 sw=4 tw=0:
