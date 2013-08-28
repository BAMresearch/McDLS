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
    lower, upper = 0., 1.

    @classmethod
    def get(cls, count = 1):
        rs = 10**(numpy.random.uniform(cls.lower, cls.upper, count))
        rs = (rs - 1) / (10**(cls.upper - cls.lower))
        return rs

class RandomExponential1(RandomExponential):
    """Alias class for RandomExponential"""
    pass

class RandomExponential2(RandomExponential):
    """Picks values with inverse logarithmic probability over )0, 1(
    , as if it were spanning two decades."""
    upper = 2.

class RandomExponential3(RandomExponential):
    """Picks values with inverse logarithmic probability over )0, 1(
    , as if it were spanning three decades."""
    upper = 3.

# vim: set ts=4 sts=4 sw=4 tw=0:
