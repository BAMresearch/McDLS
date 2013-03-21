# -*- coding: utf-8 -*-
# mcsas_test.py

import os.path
import pickle
import numpy
# requires https://nose.readthedocs.org/en/latest/
import nose

from McSAS import McSAS

TOL = 1e-16 # usual eps value for 64bit float precision
FN_TEST_DATA = "test_data.pydat"
FN_RESULT_DATA = "test_data_result.pydat"

def isEqualFloat(a, b, tol = TOL):
    """Return TRUE if both float arrays can be considered as equal.
       *tol* Tolerance in relative mean difference
       Supposed to be symmetrical."""
    # equal = numpy.allclose(a, b, rtol = 0.0, atol = TOL)

    # one could also think about using the relative mean difference:
    # http://stackoverflow.com/a/4029392
    diff = (abs(a-b) / (abs(a)+abs(b)) / 2.)
    equal = numpy.all(diff < tol)

    if not equal:
        print "value"
        print a
        print "expectation"
        print b
        print "diff"
        print diff
        print "max difference:", numpy.max(diff)
    return equal

def test():
    """Testing the algorithm in 1D. Atm, we just test as much as possible.
    Testing post-processing routines should be separated later
    as it's deterministic and though easier to test."""

    # read test data file
    with open(FN_TEST_DATA) as fd:
        qie = pickle.load(fd)

    # separate input data vectors
    q, i, e = qie[0,:], qie[1,:], qie[2,:]

    # TODO: extract settings from FN_RESULT_DATA
    # run the monte carlo routine
    mcsas = McSAS(Q = q, I = i, IERR = numpy.maximum(0.01*i, e),
                  Ncontrib = 200, Nreps = 4,
                  Convcrit = 1, Maxiter = 1e5,
                  Histscale = 'log', drhosqr = 1e30,
                  Plot = False)

    result = mcsas.result[0]
    for key, item in result.iteritems():
        print key,
        try:
            print item.shape
        except:
            pass

    expectation = None
    if not os.path.exists(FN_RESULT_DATA):
        print("Storing test results to file ... ")
        with open(FN_RESULT_DATA, 'w') as fd:
            pickle.dump(result, fd)
        print("done.")
        return
    else:
        print("Reading test result expectations from file ... ")
        with open(FN_RESULT_DATA) as fd:
            expectation = pickle.load(fd)
        print("done.")

    assert expectation is not None
    # testing against data from file
    # test only items which are averaged over all repetitions
    # individual tolerances because of large deviations for few repetitions(4)
    for key, tol in (("Hx", TOL), ("Hmid", TOL), ("Hwidth", TOL),
                     ("Hmean", 0.25), ("Hnmean", 0.25),
                     ("Hstd", 0.5), ("Hnstd", 0.5),
                     ("vfminbins", 0.1), ("nfminbins", 0.1),
                     ("Qfit", TOL), ("Imean", 0.005), ("Istd", 0.5)):
        assert isEqualFloat(result[key], expectation[key], tol),\
            "Test for {0} failed!".format(key)

if __name__ == "__main__":
    test()

# vim: set ts=4 sts=4 sw=4 tw=0:
