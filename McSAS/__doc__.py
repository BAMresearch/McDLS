# -*- coding: utf-8 -*-
# __doc__.py

r"""
Overview
========
A class and supplementary functions for Monte-Carlo fitting of SAXS patterns.
It is released under a `Creative Commons CC-BY-SA license
<http://creativecommons.org/licenses/by-sa/3.0/>`_.
Please cite as::

    Brian R. Pauw et al., J. Appl. Cryst. 46, (2013), pp. 365--371
        doi: http://dx.doi.org/10.1107/S0021889813001295

Classes and Functions Defined in This File
------------------------------------------

 - :py:class:`McSAS() <McSAS.McSAS>`:
   A class containing all the Functions required to perform a
   Monte Carlo analysis on small-angle scattering data.
    
Made possible with help from (amongst others)
---------------------------------------------

 - | Ingo Bressler <ingo.bressler@bam.de>
   | Code cleanup, modification and documentation
 - | Pawel Kwasniewski <kwasniew@esrf.fr>
   | Code cleanup and documentation
 - | Samuel Tardif
   | Derivations (mostly observability) and checking of mathematics
 - | Jan Skov Pedersen
   | checking of mathematics

A Note on Units
---------------

Internally, all length units are in meters, 
all angle units in degrees clockwise from top. 
*Intensity* is in :math:`\left[ 1 \over {m \cdot sr} \right]`,
*q* in :math:`\left[ 1 \over m \right]`.
The electron density contrast squared, *DeltaRhoSquared* is in 
    :math:`\left[ m^{-4} \right]`.
Other units may be used, but if absolute units are supplied and absolute
volume fractions required, meters are required.

Example Usage
-------------

*For detailed usage, please see the* :doc:`quickstart`
*this section will be updated shortly

Module Documentation
====================
"""

# vim: set ts=4 sts=4 sw=4 tw=0:
