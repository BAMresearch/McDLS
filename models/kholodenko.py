# -*- coding: utf-8 -*-
# models/kholodenko.py

import logging
import time
import numpy
from numpy import pi
from mcsas.McSAS import ScatteringModel
from cutesnake.algorithm import Parameter, RandomUniform
from scatteringmodel import ScatteringModel

# parameters must not be inf

from scipy.special import j1 as bessel_j1
from scipy.integrate import quad
from multiprocessing import Process, Queue, cpu_count

LASTMSG = set()

def calcParamsVsQ(func, vec, args, start = None, end = None, queue = None):
    if start is None:
        start = 0
    if end is None:
        end = len(args)
    result = numpy.empty((len(vec), end-start))
    for idx in range(start, end):
        # queue it here directly? what about non-parallel version?
        result[:, idx-start] = func(args[idx], vec)
    if queue is None:
        return result
    queue.put((start, end, result))
    queue.close()

def calcParamsVsQParallel(func, vec, args):
    if len(args) == 1:
        return calcParamsVsQ(func, vec, args)
    t0 = time.time()
    dataPerCore = max(len(args) / int(cpu_count() * 1.5), 1)
    result = numpy.zeros((len(vec), len(args)))
    queue = Queue()
    jobs = []
    for start in range(0, len(args), dataPerCore):
        end = min(start + dataPerCore, len(args))
        p = Process(target = calcParamsVsQ, args = (func, vec, args, start, end, queue))
        jobs.append(p)
    print "starting"
    for job in jobs:
        job.start()
    print "joining"
    for job in jobs:
        job.join()
    print "get results", result.shape, queue.qsize()
    resSum = 0
    while resSum < result.shape[1]:
        try:
            start, end, res = queue.get()
        except IOError:
            continue
        resSum += res.shape[1]
        print "resSum", resSum, result.shape, res.shape
        result[:, start:end] = res
    logging.info("Calculated p0 by {0} processes in {1} secs."
                 .format(len(jobs), time.time()-t0))
    del jobs[:], p, job
    return result

def core(z, qValue, kuhnLength, x):
    if z <= 0.0 or x <= 0.0:
        return 1.0
    ratio = 3.0 / kuhnLength
    if qValue < ratio:
        e = numpy.sqrt(1.0 - qValue*qValue*kuhnLength*kuhnLength/9.)
        fz = numpy.sinh(e*z) / (e*numpy.sinh(z))
    elif qValue > ratio:
        f = numpy.sqrt(qValue*qValue*kuhnLength*kuhnLength/9. - 1.0)
        fz = numpy.sin(f*z)  / (f*numpy.sinh(z))
    else: # qValue == ratio
        fz = z / numpy.sinh(z)
    res = fz * (2./x) * (1.0 - z/x)
    return res

def coreIntegral(qValue, kuhnLength, x):
    res = quad(core, 0, x,
               args = (qValue, kuhnLength, x),
               limit = 10000, full_output = 1, epsabs = 0.0, epsrel = 1e-10)
    if len(res) > 3:
        LASTMSG.add(res[-1])
    # print qValue, res[0:2]
    return res[0:2]
vectorizedCoreIntegral = numpy.vectorize(coreIntegral)

def coreIntegralOverQ(constants, qVector):
    kuhnLength, x = constants
    res = vectorizedCoreIntegral(qVector, kuhnLength, x)
    return res[0]

def calcPcs(u):
    if u <= 0.0:
        return 1.0
    res = 4. * bessel_j1(u)*bessel_j1(u) / (u*u)
    return res
vectorizedPcs = numpy.vectorize(calcPcs)

class Kholodenko(ScatteringModel):
    r"""Form factor of a worm-like structure after [Kholodenko93]_

    .. [Kholodenko93] `A. L. Kholodenko. Analytical calculation of the 
        scattering function for polymers of arbitrary flexibility using the
        dirac propagator. Macromolecules, 26:4179–4183, 1993.
        <http://dx.doi.org/10.1021/ma00068a017>`_
    """
    shortName = "Kholodenko Worm"
    parameters = (
            Parameter("radius", 1.0,
                    valueRange = (0., numpy.inf), suffix = "nm"),
            Parameter("lenKuhn", 1.0,
                    displayName = "kuhn length",
                    valueRange = (0., numpy.inf), suffix = "nm"),
            Parameter("lenContour", 1.0,
                    displayName = "contour length",
                    valueRange = (0., numpy.inf), suffix = "nm")
    )
    parameters[0].isActive = True
    parameters[1].isActive = True
    parameters[2].isActive = True

    def __init__(self):
        ScatteringModel.__init__(self)
        # some presets
        self.radius.setValueRange((1, 5))
        self.lenKuhn.setValueRange((10, 50))
        self.lenContour.setValueRange((100, 1000))

    def updateParamBounds(self, bounds):
        bounds = ScatteringModel.updateParamBounds(self, bounds)
        if len(bounds) < 1:
            return
        logging.info("bounds [pi/qmin, pi/qmax]: '{}'".format(bounds))
        return
        if len(bounds) == 1:
            logging.warning("Only one bound provided, "
                            "assuming it denotes the maximum.")
            bounds.insert(0, self.radius.valueRange[0])
        elif len(bounds) > 2:
            bounds = bounds[0:2]
            logging.info("Updating lower and upper contribution parameter bounds "
                         "to: ({0}, {1}).".format(bounds[0], bounds[1]))
            self.radius.valueRange = (min(bounds), max(bounds))

    def ff(self, dataset, paramValues):
        assert ScatteringModel.ff(self, dataset, paramValues)

        # vectorized data and arguments
        q = dataset[:, 0]
        radius = numpy.array((self.radius.value(),))
        if self.radius.isActive:
            radius = paramValues[:, 0]
        lenKuhn = numpy.array((self.lenKuhn.value(),))
        if self.lenKuhn.isActive:
            idx = int(self.radius.isActive)
            lenKuhn = paramValues[:, idx]
        lenContour = numpy.array((self.lenContour.value(),))
        if self.lenContour.isActive:
            idx = int(self.radius.isActive) + int(self.lenKuhn.isActive)
            lenContour = paramValues[:, idx]
        qr = numpy.outer(q, radius) # a matrix usually
        pcs = vectorizedPcs(qr)

        x = numpy.divide(3. * lenContour, lenKuhn)
        if len(lenKuhn) == 1:
            value = lenKuhn[0]
            lenKuhn = lenContour.copy()
            lenKuhn.fill(value)
        constants = numpy.array((lenKuhn, x)).T
        p0 = calcParamsVsQ(coreIntegralOverQ, q, constants)
        if len(LASTMSG):
            logging.warning("\n".join(["numpy.quad integration messages: "] + list(LASTMSG)))
        return p0 * pcs

    def vol(self, paramValues, compensationExponent = None):
        assert ScatteringModel.vol(self, paramValues)
        if compensationExponent is None:
            compensationExponent = self.compensationExponent
        radius = numpy.array((self.radius.value(),))
        if self.radius.isActive:
            radius = paramValues[:, 0]
        lenContour = numpy.array((self.lenContour.value(),))
        if self.lenContour.isActive:
            idx = int(self.radius.isActive) + int(self.lenKuhn.isActive)
            lenContour = paramValues[:, idx]
        volume = pi * radius*radius * lenContour
        if len(volume) <= 1 and len(volume) < len(paramValues):
            # duplicates of single value for dimension compatibility
            volume = numpy.ones(len(paramValues)) * volume
        return volume**compensationExponent

Kholodenko.factory()

if __name__ == "__main__":
    from cutesnake.datafile import PDHFile, AsciiFile
    pf = PDHFile("sasfit_kho-1-10-1000.dat")
    model = Kholodenko()
    model.radius.setValue(1.)
    model.radius.isActive = False
    model.lenKuhn.setValue(10.)
    model.lenKuhn.isActive = False
    model.lenContour.setValue(1000.)
    model.lenContour.isActive = False
    intensity = model.ff(pf.data, None).reshape(-1)
    q = pf.data[:, 0]
    oldInt = pf.data[:, 1]
    delta = abs(oldInt - intensity)
    result = numpy.dstack((q, intensity, delta))[0]
    AsciiFile.writeFile("kho.dat", result)
    # call it like this:
    # PYTHONPATH=..:../mcsas/ python brianpauwgui/kholodenko.py && gnuplot -p -e 'set logscale xy; plot "kho.dat" using 1:2:3 with errorbars'

# vim: set ts=4 sts=4 sw=4 tw=0:
