.. Find the reST syntax at http://sphinx-doc.org/rest.html
.. TODO: download mathjax in conf.py and put it into the static dir
.. http://www.sphinx-doc.org/en/stable/ext/math.html

*********************************
The Math Behind
*********************************

====
SAXS
====

The form factor is defined as:

.. math::

    \def\v#1{{\bf #1}}
    ff_{sph}(q, r) &= { 3 ~ sin(qr) - qr ~ cos(qr) \over (qr)^3 } \\
    v_{sph}(r) &= {\tiny {4 \over 3}} ~ \pi ~ r^3 \\
    v_{sph,abs}(r, \rho) &= \rho^2 ~ v_{sph}(r)

Where *q* is
    the scattering vector loaded from the data file and possibly preprocessed,
    respectively filtered by defining min/max *q* or masking invalid values
    equal or below zero.

*r* denotes
    the radius of the sphere set in the user interface (UI) or varied during
    optimization.

:math:`\rho` denotes
    the scattering length density constant of the model defined in the UI.

=============
example test
=============

.. math:: \v{e}^{i\pi} + 1 = 0
   :label: euler

Euler's identity, equation :eq:`euler`, was elected one of the most
beautiful mathematical formulas.


.. vim: set ts=4 sts=4 sw=4 tw=0:
