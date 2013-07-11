# -*- coding: utf-8 -*-
# numbergenerator.py

from abc import ABCMeta, abstractmethod
import numpy

class NumberGenerator(object):
    """Base class for number generators.
    Generates numbers in the interval [0, 1].
    Scaling is supposed to happen elsewhere."""
    __metaclass__ = ABCMeta

    @classmethod
    @abstractmethod
    def get(cls, count = 1):
        raise NotImplementedError

class RandomUniform(NumberGenerator):
    @classmethod
    def get(cls, count = 1):
        return numpy.random.uniform(size = count)

class RandomExponential(NumberGenerator):
    @classmethod
    def get(cls, count = 1):
        #numpy.random.exponential is an unbound function, may cause issues
        #as it may assume values >1. Therefore, a new sample is drawn upon its
        #occurrence
        Rs=numpy.random.exponential(size = count)
        if count>1:
            while Rs.max()>1.:
                Rs[Rs>1.]=numpy.random.exponential(size=(Rs>1.).sum())
        else:
            while Rs>1.:
                Rs=numpy.random.exponential()
        return Rs

# vim: set ts=4 sts=4 sw=4 tw=0:
