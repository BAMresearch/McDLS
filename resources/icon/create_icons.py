#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# create_icons.py

"""
Creates a set of icon PNG images at preconfigured resolutions from an input
SVG file.

From the McSAS source tree, run it like that:

    PYTHONPATH=. python3 resources/icon/create_icons.py resources/icon/mcsas.svg

"""

import sys
import re
import os
import math

sys.path.append(os.getcwd())

from gui.version import version

ICONSIZES = (16, 32, 48, 64, 128, 256, 512, 1024)

def renderSVGQt(contents, outfn, size):
    """Has font problems, letter spacing too broad, version number cropped
    (compared to inkscape output).
    Unused, but kept for reference.
    """
    from QtSvg import QSvgRenderer
    from QtCore import QByteArray
    from QtWidgets import QImage, QPainter, QApplication
    app = QApplication([])
    svg = QSvgRenderer(QByteArray(contents.encode("utf-8")))
    img = QImage(size, size, QImage.Format_ARGB32)
    painter = QPainter(img)
    svg.render(painter)
    painter.end()
    img.save(outfn)

def renderSVGrsvg(contents, outfn, size):
    """Not supported in Python3. Unused, but kept for reference."""
    import cairo
    import rsvg
    img = cairo.ImageSurface(cairo.FORMAT_ARGB32, size, size)
    ctx = cairo.Context(img)
    handle = rsvg.Handle(None, contents)
    handle.render_cairo(ctx)
    img.write_to_png(outfn)

def renderCairoSVG(contents, outfn, scale):
    """Works best with v2 and Python3. v1.0.22 for Python2 does not draw the
    gradient. Python3 is required."""
    if sys.version_info[0] < 3:
        sys.stdout.write("Sorry, Python 3.x! "
                         "With python2, gradients get not drawn properly.\n")
        sys.exit(1)
    import cairosvg
    cairosvg.svg2png(bytestring = contents.encode("utf-8"),
                     scale = scale,
                     write_to = outfn)

def adjust(contents, size):
    """Adjust the SVG code to improve output for small icon sizes."""
    # does not work due to scaling at other places in the SVG
#    contents = re.sub("(width=)\"([0-9]+)\"", "\\1\"{}\"".format(size), contents, 1)
#    contents = re.sub("(height=)\"([0-9]+)\"", "\\1\"{}\"".format(size), contents, 1)
    # get the document size, report it back for proper scaling
    width, height = None, None
    matches = re.findall("(width=)\"([0-9]+)\"", contents)
    if len(matches):
        width = int(matches[0][-1])
    matches = re.findall("(height=)\"([0-9]+)\"", contents)
    if len(matches):
        height = int(matches[0][-1])
    # update colors
    colors = version.colors()
    if len(colors) > 3:
        innerColor, outerColor, cirlesColor, shellColor = colors
        for num, color in (26, innerColor), (28, outerColor):
            contents = re.sub("(#stop{}".format(num) +
                "\\s+\\{[^\\}]*stop-color:\\s*)#[0-9a-zA-Z]+(\\s*;)",
                "\\g<1>{}\\g<2>".format(color), contents, 1)
        for elem, color in ("fill", cirlesColor), ("stroke", shellColor):
            contents = re.sub("(#shell\\s+\\{[^\\}]*"
                    + elem + ":\\s*)#[0-9a-zA-Z]+(\\s*;)",
                "\\g<1>{}\\g<2>".format(color), contents, 1)
    # set the program name
    nameFront, nameBack = version.name()[0:2], version.name()[2:]
    contents = re.sub(
        "(id=\"programNameFront\"\\s+.+>)[a-zA-Z]+(<[^>]+></text>)",
        "\\g<1>{}\\g<2>".format(nameFront), contents, 1)
    contents = re.sub(
        "(id=\"programNameBack\"\\s+.+>)[a-zA-Z]+(<[^>]+></text>)",
        "\\g<1>{}\\g<2>".format(nameBack), contents, 1)
    # set the current version number
    contents = re.sub(
        "(id=\"versionNumber\"\\s+.+>)[0-9]+\\.[0-9]+(<[^>]+></text>)",
        "\\g<1>{}\\g<2>".format(version.number()), contents, 1)
    # adjust version number size according to target icon size
    pattern = re.compile(
        "(id=\"versionNumber\"\\s+transform=\"matrix\\()([0-9\\.,-]+)(\\)\")")
    match = pattern.search(contents)
    if len(match.groups()) > 2:
        matrix = [float(v) for v in match.group(2).split(',')]
        if size < 64:
            ratio = max(width, height) / size
            scale = math.log10(ratio)
            matrix[0] *= scale
            matrix[3] *= scale
            # tweak the dx/dy adjustment for proper placement
            offset = matrix[-1] * (scale - 1.) * .2
            matrix[-2] -= offset * 1.5
            matrix[-1] -= offset
            matrix = ",".join(["{0:.3f}".format(v) for v in matrix])
            contents = pattern.sub("\\g<1>{}\\g<3>".format(matrix), contents, 1)
    if size <= 16:
        # for tiny icons hide the version number and the outer shell
        for elem in ("shell", "versionNumber"):
            contents = re.sub("(id=\"{}\")".format(elem),
                              "\\1 style=\"display:none;\"", contents, 1)
        # at the same time, enlarge the core and the label
        contents = re.sub("(id=\"core\")",
                          "\\1 transform=\"matrix(1.15,0,0,1.15)\"",
                          contents, 1)
    return contents, width, height

def createIcons(svgfn):
    fnbase = os.path.splitext(svgfn)[0]
    svgData = None
    with open(svgfn, 'rb') as fd:
        svgData = fd.read()
    contents = " ".join([l.decode("utf-8").strip() for l in svgData.splitlines()])
    for iconsize in ICONSIZES:
        outfn = "{fn}_{px:04}.png".format(fn = fnbase, px = iconsize)
        outSVG, width, height = adjust(contents, iconsize)
        renderCairoSVG(outSVG, outfn, iconsize / float(max(width, height)))

if __name__ == "__main__":
    assert len(sys.argv) > 1, (
        "Please provide an icon SVG file!")
    createIcons(sys.argv[1])

# vim: set ts=4 sw=4 sts=4 tw=0:
